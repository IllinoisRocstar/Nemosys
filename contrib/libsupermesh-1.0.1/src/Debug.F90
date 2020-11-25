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

! The following code is derived from
! debug/Debug.F90 in Fluidity git revision
! 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
!
! Fluidity copyright information (note that AUTHORS mentioned in the following
! has been renamed to Fluidity_AUTHORS):
!
!  Copyright (C) 2006 Imperial College London and others.
!  
!  Please see the AUTHORS file in the main source directory for a full list
!  of copyright holders.
!
!  Prof. C Pain
!  Applied Modelling and Computation Group
!  Department of Earth Science and Engineering
!  Imperial College London
!
!  amcgsoftware@imperial.ac.uk
!  
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation,
!  version 2.1 of the License.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!  USA

#include "libsupermesh_debug.h"

module libsupermesh_debug

  use iso_c_binding, only : c_char, c_int, c_null_char
  use iso_fortran_env, only : error_unit
  use mpi
  
  implicit none
  
  private
  
  public :: abort_pinpoint

  interface print_backtrace
    subroutine libsupermesh_print_backtrace(max_size) bind(c)
      use iso_c_binding, only : c_int
      implicit none
      integer(kind = c_int), value :: max_size
    end subroutine libsupermesh_print_backtrace
  end interface print_backtrace
  
contains

  subroutine abort_pinpoint(error, file, line_number)
    character(len = *), intent(in) :: error
    character(len = *), intent(in) :: file
    integer, intent(in) :: line_number
    
    integer :: ierr
    logical :: mpi_init
    
    call MPI_Initialized(mpi_init, ierr)
    write(error_unit, "(a)")      "*** libsupermesh error ***"
    write(error_unit, "(a,i0,a)") "Source location: (" // trim(file) // ",", line_number, ")"
    write(error_unit, "(a)")      "Error message: " // trim(error)
    call print_backtrace(max_size = 64)
    if(mpi_init) call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    stop 1
    
  end subroutine abort_pinpoint
  
  subroutine cabort_pinpoint(error, file, line_number) bind(c, name = "libsupermesh_abort_pinpoint")
    character(kind = c_char), dimension(*), intent(in) :: error
    character(kind = c_char), dimension(*), intent(in) :: file
    integer(kind = c_int), value, intent(in) :: line_number
    
    call abort_pinpoint(from_cstring(error), from_cstring(file), line_number)
    
  end subroutine cabort_pinpoint
  
  pure function len_cstring(cs) result(l)
    character(kind = c_char), dimension(*), intent(in) :: cs
    
    integer :: l
    
    l = 0
    do while(cs(l + 1) /= c_null_char)
      l = l + 1
    end do
    
  end function len_cstring
  
  pure function from_cstring(cs) result(fs)
    character(kind = c_char), dimension(*), intent(in) :: cs
    
    character(len = len_cstring(cs)) :: fs
    
    integer :: i
    
    forall(i = 1:len(fs))
      fs(i:i) = cs(i)
    end forall
    
  end function from_cstring

end module libsupermesh_debug
