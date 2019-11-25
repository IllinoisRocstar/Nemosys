!  Copyright (C) 2016-2017 The University of Edinburgh
!
!  The file is part of libsupermesh
!    
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation;
!  version 2.1 of the License.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include "libsupermesh_debug.h"

module libsupermesh_precision

#ifdef LIBSUPERMESH_DOUBLE_PRECISION
  use mpi, only : mpi_real_kind => MPI_REAL8
  use iso_fortran_env, only : real_kind => real64
#else
  use mpi, only : mpi_real_kind => MPI_REAL4
  use iso_fortran_env, only : real_kind => real32
#endif

  implicit none  

  private
  
  public :: mpi_real_kind, real_format, real_kind
  public :: compensated_sum, initialise, add, sum
  
#ifdef LIBSUPERMESH_DOUBLE_PRECISION
  character(len = *), parameter :: real_format = "es24.16e3"  
#else
  character(len = *), parameter :: real_format = "es15.8e2"
#endif

  type compensated_sum
    real(kind = real_kind) :: sum
    real(kind = real_kind) :: error
  end type compensated_sum
  
  interface initialise
    module procedure initialise_compensated_sum
  end interface initialise
  
  interface add
    module procedure add_compensated_sum
  end interface add
  
  interface sum
    module procedure sum_compensated_sum
  end interface sum
  
contains

  pure subroutine initialise_compensated_sum(cs)
    type(compensated_sum), intent(out) :: cs
    
    cs%sum = 0.0_real_kind
    cs%error = 0.0_real_kind
  
  end subroutine initialise_compensated_sum
  
  pure subroutine add_compensated_sum(cs, summand)
    type(compensated_sum), intent(inout) :: cs
    real(kind = real_kind), intent(in) :: summand
    
    real(kind = real_kind) :: compensated_summand, sum
    
    ! Kahan summation from page 9-4 of:
    !   Implementation of algorithms, Part I, Technical Report 20, W. Kahan,
    !   Lecture notes by W. S. Haugeland and D. Hough, Department of Computer
    !   Science, University of California, Berkeley, 1973
    ! See also:
    !   Pracniques: Further remarks on reducing truncation errors, W. Kahan,
    !   Communications of the ACM vol. 8, 1965, p. 40
    !
    !   The accuracy of floating point summation, N. J. Higham, SIAM Journal on
    !   Scientific Computing vol. 14, 1993, pp. 783--799
    
    compensated_summand = cs%error + summand   ! Compensated summand, including previous error
    sum = cs%sum + compensated_summand         ! New sum, including the compensated summand
    
    cs%error = cs%sum - sum                    ! Negative of what was actually added
    cs%error = cs%error + compensated_summand  ! What we failed to add
    
    cs%sum = sum
  
  end subroutine add_compensated_sum
  
  pure function sum_compensated_sum(cs) result(sum)
    type(compensated_sum), intent(in) :: cs
    
    real(kind = real_kind) :: sum
    
    ! Apply a final correction. See:
    !   Implementation of algorithms, Part I, Technical Report 20, W. Kahan,
    !   Lecture notes by W. S. Haugeland and D. Hough, Department of Computer
    !   Science, University of California, Berkeley, 1973
    !
    !   The accuracy of floating point summation, N. J. Higham, SIAM Journal on
    !   Scientific Computing vol. 14, 1993, pp. 783--799
    sum = cs%sum + cs%error
    
  end function sum_compensated_sum
  
end module libsupermesh_precision
