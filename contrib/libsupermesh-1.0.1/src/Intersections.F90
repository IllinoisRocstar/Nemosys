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

module libsupermesh_intersections

  implicit none
  
  private
  
  public :: intersections, deallocate, intersections_to_csr_sparsity
  
  type intersections
    integer, dimension(:), pointer :: v
    integer :: n
  end type intersections

  interface deallocate
    module procedure deallocate_intersections
  end interface deallocate
  
contains

  pure elemental subroutine deallocate_intersections(ints)
    type(intersections), intent(inout) :: ints

    deallocate(ints%v)

  end subroutine deallocate_intersections

  pure subroutine intersections_to_csr_sparsity(ints, indices, indptr)
    type(intersections), dimension(:), intent(in) :: ints
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: indices
    integer, dimension(:), intent(out) :: indptr

    integer :: i, j, n

    n = 0
    do i = 1, size(ints)
      n = n + ints(i)%n
    end do

    allocate(indices(n))
    indptr(1) = 1
    do i = 1, size(ints)
      indptr(i + 1) = indptr(i) + ints(i)%n
      do j = 1, ints(i)%n
        indices(indptr(i) + j - 1) = ints(i)%v(j)
      end do
    end do

  end subroutine intersections_to_csr_sparsity
  
end module libsupermesh_intersections
