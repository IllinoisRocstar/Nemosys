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

! This uses modified versions of code from femtools/Adjacency_Lists.F90 in
! Fluidity 4.1.11, first added 2015-07-05
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

module libsupermesh_graphs

  use libsupermesh_debug, only : abort_pinpoint

  implicit none

  private

  public :: eelist_type, deallocate, mesh_nelist, mesh_eelist, sloc, nneigh_ele

  type eelist_type
    integer, dimension(:, :), pointer :: v
    integer, dimension(:), pointer :: n
  end type eelist_type

  interface deallocate
    module procedure deallocate_eelist
  end interface deallocate
  
contains

  pure subroutine mesh_nelist(nnodes, enlist, nelist_indices, nelist_indptr)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: nelist_indices
    ! nnodes + 1
    integer, dimension(:), intent(out) :: nelist_indptr
    
    integer :: i, j, loc, nelements
    integer, dimension(:), allocatable :: nelist_n

    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    ! Construct the node-element list. Based on NODELE in Fluidity
    ! femtools/Adjacency_Lists.F90 (see e.g. Fluidity version 4.1.11). Code
    ! first added 05/07/2015.
    allocate(nelist_n(nnodes))
    nelist_n = 0
    do i = 1, nelements
      do j = 1, loc
        nelist_n(enlist(j, i)) = nelist_n(enlist(j, i)) + 1
      end do
    end do
    nelist_indptr(1) = 1
    do i = 1, nnodes
      nelist_indptr(i + 1) = nelist_indptr(i) + nelist_n(i)
    end do
    allocate(nelist_indices(nelist_indptr(nnodes + 1) - 1))
    nelist_n = 0
    do i = 1, nelements
      do j = 1, loc
        nelist_indices(nelist_indptr(enlist(j, i)) + nelist_n(enlist(j, i))) = i
        nelist_n(enlist(j, i)) = nelist_n(enlist(j, i)) + 1
      end do
    end do
    deallocate(nelist_n)
  
  end subroutine mesh_nelist    

  subroutine mesh_eelist(nnodes, enlist, sloc, eelist, nneigh_ele)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    integer, intent(in) :: sloc
    type(eelist_type), intent(out) :: eelist
    integer, optional, intent(in) :: nneigh_ele

    integer :: ele, i, j, k, k_max, k_min, loc, nelements, lnneigh_ele

    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable :: nelist_indices, nelist_indptr

    integer, dimension(:), allocatable :: seen, seen_eles
    integer :: nseen_eles

    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    allocate(nelist_indptr(nnodes + 1))
    call mesh_nelist(nnodes, enlist, nelist_indices, nelist_indptr)
    
    allocate(seen_eles(nelements), seen(nelements))
    seen = 0
    if(present(nneigh_ele)) then
      lnneigh_ele = nneigh_ele
    else
      lnneigh_ele = loc
    end if
    allocate(eelist%v(lnneigh_ele, nelements), eelist%n(nelements))
    eelist%n = 0
    do i = 1, nelements
      nseen_eles = 0
      loc_loop: do j = 1, loc
        k_min = nelist_indptr(enlist(j, i))
        k_max = nelist_indptr(enlist(j, i) + 1) - 1
        do k = k_min, k_max
          ele = nelist_indices(k)
          if(ele /= i) then
            seen(ele) = seen(ele) + 1
            if(seen(ele) == 1) then
              nseen_eles = nseen_eles + 1
              seen_eles(nseen_eles) = ele
            end if
            if(seen(ele) == sloc) then
              eelist%n(i) = eelist%n(i) + 1
#ifndef NDEBUG
              if(eelist%n(i) > lnneigh_ele) then
                libsupermesh_abort("Invalid connectivity")
              end if
#endif
              eelist%v(eelist%n(i), i) = ele
#ifdef NDEBUG
              if(eelist%n(i) == lnneigh_ele) exit loc_loop
#else
            else if(seen(ele) > sloc) then
              libsupermesh_abort("Invalid connectivity")
#endif
            end if
          end if
        end do
      end do loc_loop
      do j = 1, nseen_eles
        seen(seen_eles(j)) = 0
      end do
    end do
    deallocate(nelist_indices, nelist_indptr, seen_eles, seen)

  end subroutine mesh_eelist

  pure elemental function sloc(dim, loc)
    integer, intent(in) :: dim
    integer, intent(in) :: loc

    integer :: sloc

    if(loc == dim + 1) then
      ! Simplex
      sloc = loc - 1
    else if(loc == 2 ** dim) then
      ! Cubical
      sloc = loc / 2
    else
      sloc = -1
    end if

  end function sloc
  
  pure elemental function nneigh_ele(dim, loc)
    integer, intent(in) :: dim
    integer, intent(in) :: loc

    integer :: nneigh_ele

    if(loc == dim + 1) then
      ! Simplex
      nneigh_ele = loc
    else if(loc == 2 ** dim) then
      ! Cubical
      nneigh_ele = 2 * dim
    else
      nneigh_ele = -1
    end if
  
  end function nneigh_ele

  pure elemental subroutine deallocate_eelist(eelist)
    type(eelist_type), intent(inout) :: eelist

    deallocate(eelist%v, eelist%n)

  end subroutine deallocate_eelist

end module libsupermesh_graphs
