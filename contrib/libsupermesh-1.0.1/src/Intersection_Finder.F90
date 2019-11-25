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
! femtools/Intersection_finder.F90 in Fluidity git revision
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

module libsupermesh_intersection_finder

  use iso_c_binding, only : c_double, c_int, c_ptr
!   use iso_fortran_env, only : error_unit
  
  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_graphs, only : eelist_type, deallocate, mesh_eelist, sloc, &
    & nneigh_ele
  use libsupermesh_intersections, only : intersections, deallocate, &
    & intersections_to_csr_sparsity
  use libsupermesh_octree_intersection_finder, only : octree_node, &
    & octree_type, deallocate, octree_intersection_finder, allocate, query
  use libsupermesh_precision, only : real_kind
  use libsupermesh_quadtree_intersection_finder, only : quadtree_node, &
    & quadtree_type, deallocate, quadtree_intersection_finder, allocate, query

  implicit none

  private
  
  interface cbuild_rtree
    subroutine libsupermesh_build_rtree(rtree, dim, nnodes, positions, loc, nelements, enlist) bind(c)
      use iso_c_binding, only : c_double, c_int, c_ptr
      implicit none
      type(c_ptr) :: rtree
      integer(kind = c_int), value :: dim
      integer(kind = c_int), value :: nnodes
      real(kind = c_double), dimension(dim, nnodes) :: positions
      integer(kind = c_int), value :: loc
      integer(kind = c_int), value :: nelements
      integer(kind = c_int), dimension(loc, nelements) :: enlist
    end subroutine libsupermesh_build_rtree
  end interface cbuild_rtree
    
  interface cquery_rtree
    subroutine libsupermesh_query_rtree(rtree, dim, loc_a, element_a, neles_b) bind(c)
      use iso_c_binding, only : c_double, c_int, c_ptr
      implicit none
      type(c_ptr) :: rtree
      integer(kind = c_int), value :: dim
      integer(kind = c_int), value :: loc_a
      real(kind = c_double), dimension(dim, loc_a) :: element_a
      integer(kind = c_int) :: neles_b
    end subroutine libsupermesh_query_rtree
  end interface cquery_rtree
  
  interface cquery_rtree_intersections
    subroutine libsupermesh_query_rtree_intersections(rtree, eles_b) bind(c)
      use iso_c_binding, only : c_int, c_ptr
      type(c_ptr) :: rtree
      integer(kind = c_int), dimension(*) :: eles_b
    end subroutine libsupermesh_query_rtree_intersections
  end interface cquery_rtree_intersections
  
  interface cdeallocate_rtree
    subroutine libsupermesh_deallocate_rtree(rtree) bind(c)
      use iso_c_binding, only : c_ptr
      type(c_ptr) :: rtree
    end subroutine libsupermesh_deallocate_rtree
  end interface cdeallocate_rtree

  public :: intersections, deallocate, intersection_finder, &
    & advancing_front_intersection_finder, rtree_intersection_finder, &
    & quadtree_intersection_finder, octree_intersection_finder, &
    & tree_intersection_finder, sort_intersection_finder, &
    & brute_force_intersection_finder
  public :: octree_node, octree_type, allocate, query
  public :: quadtree_node, quadtree_type
  public :: tree_type
  public :: rtree_type
  public :: rtree_intersection_finder_set_input, &
    & rtree_intersection_finder_find, rtree_intersection_finder_query_output, &
    & rtree_intersection_finder_get_output, rtree_intersection_finder_reset
  public :: tree_intersection_finder_set_input, &
    & tree_intersection_finder_find, tree_intersection_finder_query_output, &
    & tree_intersection_finder_get_output, tree_intersection_finder_reset

  interface intersection_finder
    module procedure intersection_finder_intersections, &
      & intersection_finder_csr_sparsity
  end interface intersection_finder

  interface advancing_front_intersection_finder
    module procedure advancing_front_intersection_finder_intersections, &
      & advancing_front_intersection_finder_csr_sparsity
  end interface advancing_front_intersection_finder

  interface rtree_intersection_finder
    module procedure rtree_intersection_finder_intersections, &
      & rtree_intersection_finder_csr_sparsity
  end interface rtree_intersection_finder

  interface tree_intersection_finder
    module procedure tree_intersection_finder_intersections, &
      & tree_intersection_finder_csr_sparsity
  end interface tree_intersection_finder
  
  interface sort_intersection_finder
    module procedure sort_intersection_finder_rank_1_intersections, &
      & sort_intersection_finder_rank_2_intersections, &
      & sort_intersection_finder_rank_1_csr_sparsity, &
      & sort_intersection_finder_rank_2_csr_sparsity
  end interface sort_intersection_finder

  interface brute_force_intersection_finder
    module procedure brute_force_intersection_finder_intersections, &
      & brute_force_intersection_finder_csr_sparsity
  end interface brute_force_intersection_finder
  
  type rtree_type
    type(c_ptr) :: rtree
  end type rtree_type
  
  type tree_type
    integer :: dim
    type(quadtree_type), pointer :: quadtree
    type(octree_type), pointer :: octree
  end type tree_type
  
  interface allocate
    module procedure allocate_tree, allocate_rtree
  end interface allocate   
  
  interface query
    module procedure query_tree_allocatable, query_tree_pointer, &
      & query_rtree_allocatable, query_rtree_pointer
  end interface query
  
  interface deallocate
    module procedure deallocate_tree, deallocate_rtree
  end interface deallocate    
    
  logical, save :: rtree_allocated = .false.
  integer, dimension(:), allocatable, save :: rtree_eles
  type(rtree_type), save :: rtree
  
  logical, save :: tree_allocated = .false.
  integer, dimension(:), allocatable, save :: tree_eles
  type(tree_type), save :: tree

contains

  pure function bbox(coords)
    ! dim x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: coords

    real(kind = real_kind), dimension(2, size(coords, 1)) :: bbox

    integer :: i, j

    bbox(1, :) = coords(:, 1)
    bbox(2, :) = coords(:, 1)
    do i = 2, size(coords, 2)
      do j = 1, size(coords, 1)
        bbox(1, j) = min(bbox(1, j), coords(j, i))
        bbox(2, j) = max(bbox(2, j), coords(j, i))
      end do
    end do

  end function bbox

  subroutine intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab
    
    select case(size(positions_a, 1))
      case(1)
        call sort_intersection_finder(positions_a(1, :), enlist_a, positions_b(1, :), enlist_b, map_ab)
      case(2:3)
        call rtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
      case default
        libsupermesh_abort("Invalid dimension")
    end select
  
  end subroutine intersection_finder_intersections
  
  subroutine intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr

    select case(size(positions_a, 1))
      case(1)
        call sort_intersection_finder(positions_a(1, :), enlist_a, positions_b(1, :), enlist_b, map_ab_indices, map_ab_indptr)
      case(2:3)
        call rtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
      case default
        libsupermesh_abort("Invalid dimension")
    end select
  
  end subroutine intersection_finder_csr_sparsity

  ! Advancing front intersection finder, as described in P. E. Farrell and
  ! J. R. Maddison, "Conservative interpolation between volume meshes by local
  ! Galerkin projection", Computer Methods in Applied Mechanics and Engineering,
  ! 200, pp. 89--100, 2011
  subroutine advancing_front_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab

    integer :: i, j
    integer :: dim, nnodes_b, nnodes_a, nelements_b, nelements_a

    real(kind = real_kind), dimension(2, size(positions_a, 1)) :: bbox_a
    real(kind = real_kind), dimension(:, :, :), allocatable :: bboxes_b
    type(eelist_type) :: eelist_b, eelist_a

    integer :: clue_a, ele_b, ele_a, loc_b, loc_a, neigh_b, neigh_a, seed_a
    integer, dimension(:), allocatable :: front_b, ints
    logical, dimension(:), allocatable :: seen_b, seen_a
    integer, dimension(:, :), allocatable :: front_a
    integer :: nfront_b, nfront_a, nints
!     integer :: nsub_a

    dim = size(positions_b, 1)
    nnodes_b = size(positions_b, 2)
    nnodes_a = size(positions_a, 2)
    loc_b = size(enlist_b, 1)
    nelements_b = size(enlist_b, 2)
    loc_a = size(enlist_a, 1)
    nelements_a = size(enlist_a, 2)
    if(nelements_b == 0 .or. nelements_a == 0) then
      do ele_a = 1, nelements_a
        allocate(map_ab(ele_a)%v(0))
        map_ab(ele_a)%n = 0
      end do
      return
    end if

    call mesh_eelist(nnodes_b, enlist_b, sloc(dim, loc_b), eelist_b, &
      & nneigh_ele = nneigh_ele(dim, loc_b))
    call mesh_eelist(nnodes_a, enlist_a, sloc(dim, loc_a), eelist_a, &
      & nneigh_ele = nneigh_ele(dim, loc_a))
    allocate(bboxes_b(2, dim, nelements_b))
    do i = 1, nelements_b
      bboxes_b(:, :, i) = bbox(positions_b(:, enlist_b(:, i)))
    end do

    allocate(seen_b(nelements_b), seen_a(nelements_a), front_b(nelements_b), &
      & front_a(2, nelements_a), ints(nelements_b))
    seen_b = .false.
    seen_a = .false.
    ! Stage 0: Initial target mesh seed element
    seed_a = 1
!     nsub_a = 0
    seed_a_loop: do
!       nsub_a = nsub_a + 1
      seen_a(seed_a) = .true.
      bbox_a = bbox(positions_a(:, enlist_a(:, seed_a)))

      ! Stage 1a: Find intersections with the target mesh seed element via a
      ! brute force search
      nints = 0
      do ele_b = 1, nelements_b
        if(bboxes_intersect(bbox_a, bboxes_b(:, :, ele_b))) then
          nints = nints + 1
          ints(nints) = ele_b
        end if
      end do
      allocate(map_ab(seed_a)%v(nints))
      map_ab(seed_a)%v = ints(:nints)
      map_ab(seed_a)%n = nints

      ! Stage 1b: Advance the target mesh front
      nfront_a = 0
      do i = 1, eelist_a%n(seed_a)
        neigh_a = eelist_a%v(i, seed_a)
        if(.not. seen_a(neigh_a)) then
          nfront_a = nfront_a + 1
          front_a(1, nfront_a) = neigh_a
          front_a(2, nfront_a) = seed_a
          seen_a(neigh_a) = .true.
        end if
      end do

      do while(nfront_a > 0)
        ele_a = front_a(1, nfront_a)
        clue_a = front_a(2, nfront_a)
        nfront_a = nfront_a - 1
        bbox_a = bbox(positions_a(:, enlist_a(:, ele_a)))

        ! Stage 2a: Initialise the donor mesh front
        nfront_b = map_ab(clue_a)%n
        front_b(:nfront_b) = map_ab(clue_a)%v
        seen_b(front_b(:nfront_b)) = .true.

        ! Stage 2b: Find intersections with the target mesh element by
        ! advancing the donor mesh front
        nints = 0
        i = 1
        do while(i <= nfront_b)
          ele_b = front_b(i)
          if(bboxes_intersect(bbox_a, bboxes_b(:, :, ele_b))) then
            ! An intersection has been found
            nints = nints + 1
            ints(nints) = ele_b
            ! Advance the donor mesh front
            do j = 1, eelist_b%n(ele_b)
              neigh_b = eelist_b%v(j, ele_b)
              if(.not. seen_b(neigh_b)) then
                nfront_b = nfront_b + 1
                front_b(nfront_b) = neigh_b
                seen_b(neigh_b) = .true.
              end if
            end do
          end if
          i = i + 1
        end do
        do i = 1, nfront_b
          seen_b(front_b(i)) = .false.
        end do
!         if(nints == 0) then
!           write(error_unit, "(a,i0)") "WARNING: Failed to find intersections for target element ", ele_a
!         end if
        allocate(map_ab(ele_a)%v(nints))
        map_ab(ele_a)%v = ints(:nints)
        map_ab(ele_a)%n = nints

        ! Stage 2c: Advance the target mesh front
        do i = 1, eelist_a%n(ele_a)
          neigh_a = eelist_a%v(i, ele_a)
          if(.not. seen_a(neigh_a)) then
            nfront_a = nfront_a + 1
            front_a(1, nfront_a) = neigh_a
            front_a(2, nfront_a) = ele_a
            seen_a(neigh_a) = .true.
          end if
        end do
      end do

      ! Stage 3: Find a new target mesh seed
      do while(seen_a(seed_a))
        seed_a = seed_a + 1
        if(seed_a > nelements_a) exit seed_a_loop
      end do
!       if(nsub_a == 1) then
!         write(error_unit, "(a)") "WARNING: Target mesh is not connected"
!       end if
    end do seed_a_loop
!     if(nsub_a > 1) then
!       write(error_unit ,"(a,i0)") "WARNING: Number of target connected sub-domains = ", nsub_a
!     end if

    deallocate(seen_b, seen_a, front_b, front_a, ints)
    call deallocate(eelist_b)
    call deallocate(eelist_a)
    deallocate(bboxes_b)
  
  contains

    pure function bboxes_intersect(bbox_1, bbox_2) result(intersect)
      ! 2 x dim
      real(kind = real_kind), dimension(:, :), intent(in) :: bbox_1
      ! 2 x dim
      real(kind = real_kind), dimension(:, :), intent(in) :: bbox_2

      logical :: intersect

      integer :: i

      do i = 1, size(bbox_1, 2)
        ! Strict inequalities required here for the advancing front intersection
        ! finder to work with axis aligned elements
        if(bbox_2(2, i) < bbox_1(1, i) .or. bbox_2(1, i) > bbox_1(2, i)) then
          intersect = .false.
          return
        end if
      end do
      intersect = .true.

    end function bboxes_intersect

  end subroutine advancing_front_intersection_finder_intersections

  subroutine advancing_front_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr

    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)

  end subroutine advancing_front_intersection_finder_csr_sparsity
  
  subroutine allocate_rtree(rtree, positions, enlist)
    type(rtree_type), intent(out) :: rtree
    ! dim x nnodes
    real(kind = real_kind), dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    
    integer(kind = c_int) :: dim, loc, nelements, nnodes
    
    dim = size(positions, 1)
    if(.not. any(dim == (/2, 3/))) then
      libsupermesh_abort("Invalid dimension")
    end if    
    nnodes = size(positions, 2)
    loc = size(enlist, 1)
    nelements = size(enlist, 2)
    
    call cbuild_rtree(rtree%rtree, dim, nnodes, real(positions, kind = c_double), loc, nelements, enlist)
  
  end subroutine allocate_rtree
  
  subroutine query_rtree_allocatable(rtree, element_a, eles_b)
    type(rtree_type), intent(inout) :: rtree
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    integer, dimension(:), allocatable, intent(out) :: eles_b
    
    integer(kind = c_int) :: dim, loc_a, neles_b
    
    dim = size(element_a, 1)
    loc_a = size(element_a, 2)    
    
    call cquery_rtree(rtree%rtree, dim, loc_a, real(element_a, kind = c_double), neles_b)
    allocate(eles_b(neles_b))
    call cquery_rtree_intersections(rtree%rtree, eles_b)
  
  end subroutine query_rtree_allocatable
    
  subroutine query_rtree_pointer(rtree, element_a, eles_b)
    type(rtree_type), intent(inout) :: rtree
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    integer, dimension(:), pointer, intent(out) :: eles_b
    
    integer(kind = c_int) :: dim, loc_a, neles_b
    
    dim = size(element_a, 1)
    loc_a = size(element_a, 2)    
    
    call cquery_rtree(rtree%rtree, dim, loc_a, real(element_a, kind = c_double), neles_b)
    allocate(eles_b(neles_b))
    call cquery_rtree_intersections(rtree%rtree, eles_b)
  
  end subroutine query_rtree_pointer
  
  subroutine deallocate_rtree(rtree)
    type(rtree_type), intent(inout) :: rtree
    
    call cdeallocate_rtree(rtree%rtree)
    
  end subroutine deallocate_rtree

  subroutine rtree_intersection_finder_set_input(positions, enlist)
    ! dim x nnodes
    real(kind = real_kind), dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist

    call rtree_intersection_finder_reset()
    call allocate(rtree, positions, enlist)
    rtree_allocated = .true.

  end subroutine rtree_intersection_finder_set_input

  subroutine rtree_intersection_finder_find(element_a)
    ! dim x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a

    assert(rtree_allocated)
    if(allocated(rtree_eles)) deallocate(rtree_eles)
    call query(rtree, element_a, rtree_eles)

  end subroutine rtree_intersection_finder_find
  
  subroutine rtree_intersection_finder_query_output(neles_b)
    integer, intent(out) :: neles_b
    
    assert(allocated(rtree_eles))
    neles_b = size(rtree_eles)
    
  end subroutine rtree_intersection_finder_query_output
  
  subroutine rtree_intersection_finder_get_output(ele_b, i)
    integer, intent(out) :: ele_b
    integer, intent(in) :: i
    
    assert(allocated(rtree_eles))
    ele_b = rtree_eles(i)
    
  end subroutine rtree_intersection_finder_get_output
  
  subroutine rtree_intersection_finder_reset()  
    if(.not. rtree_allocated) return
  
    if(allocated(rtree_eles)) deallocate(rtree_eles)
    call deallocate(rtree)
    rtree_allocated = .false.
    
  end subroutine rtree_intersection_finder_reset

  subroutine rtree_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab

    integer(kind = c_int) :: dim, ele_a, loc_a, nelements_a
    type(rtree_type) :: rtree
    
    dim = size(positions_a, 1)
    loc_a = size(enlist_a, 1)
    nelements_a = size(enlist_a, 2)
    
    call allocate(rtree, positions_b, enlist_b)
    do ele_a = 1, nelements_a
      call cquery_rtree(rtree%rtree, dim, loc_a, real(positions_a(:, enlist_a(:, ele_a)), kind = c_double), map_ab(ele_a)%n)
      allocate(map_ab(ele_a)%v(map_ab(ele_a)%n))
      call cquery_rtree_intersections(rtree%rtree, map_ab(ele_a)%v)   
    end do
    call deallocate(rtree)

  end subroutine rtree_intersection_finder_intersections

  subroutine rtree_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr

    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call rtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)

  end subroutine rtree_intersection_finder_csr_sparsity
  
  subroutine allocate_tree(tree, positions, enlist)
    type(tree_type), intent(out) :: tree
    ! dim x nnodes
    real(kind = real_kind), dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist

    tree%dim = size(positions, 1)
    select case(tree%dim)
      case(2)
        allocate(tree%quadtree)
        call allocate(tree%quadtree, positions, enlist)
      case(3)
        allocate(tree%octree)
        call allocate(tree%octree, positions, enlist)
      case default
        libsupermesh_abort("Invalid dimension")
    end select
  
  end subroutine allocate_tree
  
  subroutine query_tree_allocatable(tree, element_a, eles_b)
    type(tree_type), intent(inout) :: tree
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    integer, dimension(:), allocatable, intent(out) :: eles_b
        
    select case(tree%dim)
      case(2)
        call query(tree%quadtree, element_a, eles_b)
      case(3)
        call query(tree%octree, element_a, eles_b)
      case default
        libsupermesh_abort("Invalid dimension")
    end select
  
  end subroutine query_tree_allocatable
  
  subroutine query_tree_pointer(tree, element_a, eles_b)
    type(tree_type), intent(inout) :: tree
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    integer, dimension(:), pointer, intent(out) :: eles_b
        
    select case(tree%dim)
      case(2)
        call query(tree%quadtree, element_a, eles_b)
      case(3)
        call query(tree%octree, element_a, eles_b)
      case default
        libsupermesh_abort("Invalid dimension")
    end select
  
  end subroutine query_tree_pointer
  
  subroutine deallocate_tree(tree)
    type(tree_type), intent(inout) :: tree
    
    select case(tree%dim)
      case(2)
        call deallocate(tree%quadtree)
      case(3)
        call deallocate(tree%octree)
      case default
        libsupermesh_abort("Invalid dimension")
    end select
    
  end subroutine deallocate_tree

  subroutine tree_intersection_finder_set_input(positions, enlist)
    ! dim x nnodes
    real(kind = real_kind), dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    
    call tree_intersection_finder_reset()
    call allocate(tree, positions, enlist)
    tree_allocated = .true.
    
  end subroutine tree_intersection_finder_set_input

  subroutine tree_intersection_finder_find(element_a)
    ! dim x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    
    assert(tree_allocated)
    if(allocated(tree_eles)) deallocate(tree_eles)
    call query(tree, element_a, tree_eles)
    
  end subroutine tree_intersection_finder_find
  
  subroutine tree_intersection_finder_query_output(neles_b)
    integer, intent(out) :: neles_b
    
    assert(allocated(tree_eles))
    neles_b = size(tree_eles)
  
  end subroutine tree_intersection_finder_query_output
  
  subroutine tree_intersection_finder_get_output(ele_b, i)
    integer, intent(out) :: ele_b
    integer, intent(in) :: i
        
    assert(allocated(tree_eles))
    ele_b = tree_eles(i)
      
  end subroutine tree_intersection_finder_get_output
  
  subroutine tree_intersection_finder_reset()    
    if(.not. tree_allocated) return
  
    if(allocated(tree_eles)) deallocate(tree_eles)
    call deallocate(tree)
    tree_allocated = .false.
  
  end subroutine tree_intersection_finder_reset

  subroutine tree_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab

    select case(size(positions_a, 1))
      case(2)
        call quadtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
      case(3)
        call octree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
      case default
        libsupermesh_abort("Invalid dimension")
    end select

  end subroutine tree_intersection_finder_intersections

  subroutine tree_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr

    select case(size(positions_a, 1))
      case(2)
        call quadtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
      case(3)
        call octree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
      case default
        libsupermesh_abort("Invalid dimension")
    end select

  end subroutine tree_intersection_finder_csr_sparsity
  
  pure subroutine sort_intersection_finder_rank_1_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! nnodes_a
    real(kind = real_kind), dimension(:), intent(in) :: positions_a
    ! 2 x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! nnodes_b
    real(kind = real_kind), dimension(:), intent(in) :: positions_b
    ! 2 x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab
    
    integer :: ele_a, ele_b, i, j, nelements_a, nelements_b, nints
    integer, dimension(:), allocatable :: indices_a, indices_b, ints, work
    real(kind = real_kind), dimension(2) :: interval_a, interval_b
    real(kind = real_kind), dimension(:), allocatable :: left_positions_a, &
      & left_positions_b
    
    nelements_a = size(enlist_a, 2)
    nelements_b = size(enlist_b, 2)
    
    allocate(left_positions_a(nelements_a), left_positions_b(nelements_b))
    do ele_a = 1, nelements_a
      left_positions_a(ele_a) = minval(positions_a(enlist_a(:, ele_a)))
    end do
    do ele_b = 1, nelements_b
      left_positions_b(ele_b) = minval(positions_b(enlist_b(:, ele_b)))
    end do
    
    allocate(indices_a(nelements_a), indices_b(nelements_b), work(max(nelements_a, nelements_b)))
    call merge_sort(left_positions_a, indices_a, work(:nelements_a))
    call merge_sort(left_positions_b, indices_b, work(:nelements_b))
    deallocate(left_positions_a, left_positions_b, work)
    
    allocate(ints(nelements_b))
    j = 1
    do i = 1, nelements_a
      ele_a = indices_a(i)
      interval_a = positions_a(enlist_a(:, ele_a))
      nints = 0
      do while(j <= nelements_b)
        ele_b = indices_b(j)
        interval_b = positions_b(enlist_b(:, ele_b))
        if(maxval(interval_b) <= minval(interval_a)) then
          j = j + 1
        else if(minval(interval_b) >= maxval(interval_a)) then
          exit
        else
          nints = nints + 1
          ints(nints) = ele_b
          j = j + 1
        end if
      end do
      if(nints > 0) j = j - 1
      allocate(map_ab(ele_a)%v(nints))
      map_ab(ele_a)%v = ints(:nints)
      map_ab(ele_a)%n = nints
    end do
    
    deallocate(indices_a, indices_b, ints)
  
  contains
  
    ! A very basic merge sort implementation
    pure recursive subroutine merge_sort(v, indices, work)
      real(kind = real_kind), dimension(:), intent(in) :: v
      integer, dimension(:), intent(inout) :: indices
      integer, dimension(:), intent(out) :: work
      
      integer :: i, i_1, i_2, j, n 
      
      n = size(v, 1)
      
      if(n <= 4) then
        ! Switch to a basic quadratic sort for small inputs
        do i = 1, n
          indices(i) = i
          do j = i + 1, n
            ! < here for a stable sort
            if(v(j) < v(indices(i))) indices(i) = j
          end do
        end do
      else
        ! Otherwise split, recurse, and merge
        call merge_sort(v(1:n / 2), indices(1:n / 2), work(1:n / 2))
        call merge_sort(v((n / 2) + 1:n), indices((n / 2) + 1:n), work((n / 2) + 1:n))
        work(1:n / 2) = indices(1:n / 2)
        work((n / 2) + 1:n) = indices((n / 2) + 1:n) + (n / 2)
        i_1 = 1
        i_2 = (n / 2) + 1
        do i = 1, n
          if(i_1 > (n / 2)) then
            indices(i:n) = work(i_2:n)
            exit
          else if(i_2 > n) then
            indices(i:n) = work(i_1:n / 2)
            exit
          ! <= here for a stable sort
          else if(v(work(i_1)) <= v(work(i_2))) then
            indices(i) = work(i_1)
            i_1 = i_1 + 1
          else
            indices(i) = work(i_2)
            i_2 = i_2 + 1
          end if
        end do
      end if
      
    end subroutine merge_sort
    
  end subroutine sort_intersection_finder_rank_1_intersections
  
  pure subroutine sort_intersection_finder_rank_2_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! 1 x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! 2 x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! 1 x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! 2 x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab
    
    call sort_intersection_finder(positions_a(1, :), enlist_a, positions_b(1, :), enlist_b, map_ab)
  
  end subroutine sort_intersection_finder_rank_2_intersections

  pure subroutine sort_intersection_finder_rank_1_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! nnodes_a
    real(kind = real_kind), dimension(:), intent(in) :: positions_a
    ! 2 x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! nnodes_b
    real(kind = real_kind), dimension(:), intent(in) :: positions_b
    ! 2 x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr
    
    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call sort_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)
  
  end subroutine sort_intersection_finder_rank_1_csr_sparsity

  pure subroutine sort_intersection_finder_rank_2_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! 1 x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! 2 x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! 1 x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! 2 x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr
    
    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call sort_intersection_finder(positions_a(1, :), enlist_a, positions_b(1, :), enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)
  
  end subroutine sort_intersection_finder_rank_2_csr_sparsity

  ! Brute force intersection finder.
  pure subroutine brute_force_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab

    integer :: ele_a, ele_b, nelements_a, nelements_b
    real(kind = real_kind), dimension(2, size(positions_a, 1)) :: bbox_a, &
      & bbox_b

    integer, dimension(:), allocatable :: ints
    integer :: nints

    nelements_a = size(enlist_a, 2)
    nelements_b = size(enlist_b, 2)

    allocate(ints(nelements_b))
    do ele_a = 1, nelements_a
      bbox_a = bbox(positions_a(:, enlist_a(:, ele_a)))
      nints = 0
      do ele_b = 1, nelements_b
        bbox_b = bbox(positions_b(:, enlist_b(:, ele_b)))
        if(bboxes_intersect(bbox_a, bbox_b)) then
          nints = nints + 1
          ints(nints) = ele_b
        end if
      end do

      map_ab(ele_a)%n = nints
      allocate(map_ab(ele_a)%v(nints))
      map_ab(ele_a)%v = ints(:nints)
    end do
    deallocate(ints)
    
  contains

    pure function bboxes_intersect(bbox_1, bbox_2) result(intersect)
      ! 2 x dim
      real(kind = real_kind), dimension(:, :), intent(in) :: bbox_1
      ! 2 x dim
      real(kind = real_kind), dimension(:, :), intent(in) :: bbox_2

      logical :: intersect

      integer :: i

      do i = 1, size(bbox_1, 2)
        if(bbox_2(2, i) <= bbox_1(1, i) .or. bbox_2(1, i) >= bbox_1(2, i)) then
          intersect = .false.
          return
        end if
      end do
      intersect = .true.

    end function bboxes_intersect

  end subroutine brute_force_intersection_finder_intersections

  pure subroutine brute_force_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
    ! dim x nnodes_a
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real(kind = real_kind), dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr

    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call brute_force_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)

  end subroutine brute_force_intersection_finder_csr_sparsity

end module libsupermesh_intersection_finder
