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

module libsupermesh_halo_ownership

  use mpi

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_read_halos, only : halo_type

  implicit none

  private

  public :: element_ownership, node_ownership

contains

  subroutine element_ownership(nnodes, enlist, halo, ele_owner)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    ! level 2 halo
    type(halo_type), intent(in) :: halo
    ! nelements
    integer, dimension(:), intent(out) :: ele_owner

    integer :: i
    integer, dimension(:), allocatable :: node_owner

    allocate(node_owner(nnodes))
    call node_ownership(nnodes, halo, node_owner)
    
    do i = 1, size(enlist, 2)
      ele_owner(i) = minval(node_owner(enlist(:, i)))
    end do

    deallocate(node_owner)

  end subroutine element_ownership

  subroutine node_ownership(nnodes, halo, node_owner)
    integer, intent(in) :: nnodes
    ! level 2 halo
    type(halo_type), intent(in) :: halo
    ! nelements
    integer, dimension(:), intent(out) :: node_owner

    integer :: i, ierr, process

    call MPI_Comm_rank(halo%comm, process, ierr)
    assert(ierr == MPI_SUCCESS)
    node_owner(:halo%npnodes) = process
    if(size(node_owner) > halo%npnodes) node_owner(halo%npnodes + 1:) = huge(node_owner(halo%npnodes + 1))

    do i = 1, halo%nprocs
      if(size(halo%recv(i)%v) > 0) then
        assert(all(halo%recv(i)%v > halo%npnodes))
        assert(all(halo%recv(i)%v <= nnodes))
        node_owner(halo%recv(i)%v) = i - 1
      end if
    end do

  end subroutine node_ownership

end module libsupermesh_halo_ownership
