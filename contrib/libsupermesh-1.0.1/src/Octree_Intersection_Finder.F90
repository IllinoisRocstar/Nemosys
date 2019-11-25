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

#define TREE_DIM 3
#define TREE_NCHILDREN 8

module libsupermesh_octree_intersection_finder

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_intersections, only : intersections, deallocate, &
    & intersections_to_csr_sparsity
  use libsupermesh_precision, only : real_kind
  use libsupermesh_quadtree_intersection_finder, only : max_nelist_degree

  implicit none

  private
  
  public :: octree_type, octree_node, octree_element, max_nelist_degree, &
    & octree_intersection_finder, allocate, query, deallocate
  
  ! Octree node
  type octree_node
    ! Number of stored elements. If this is <= max_n then this is a leaf node.
    integer :: n = 0
    ! Maximum number of stored elements
    integer :: max_n  
    ! Node bounding box
    real(kind = real_kind), dimension(2, TREE_DIM) :: bbox
    ! If this is a leaf node, elements stored by this node
    type(octree_element), dimension(:), pointer :: elements
    ! Otherwise, children of this node
    type(octree_node), dimension(:), pointer :: children => null()
  end type octree_node
  
  ! Octree stored elements
  type octree_element
    ! Element bounding box
    real(kind = real_kind), dimension(2, TREE_DIM) :: bbox
    ! Element index
    integer :: ele
  end type octree_element
  
  type octree_type
    type(octree_node), pointer :: octree
    integer, dimension(:), pointer :: eles_b
    integer, pointer :: neles_b
    logical, dimension(:), pointer :: seen_ele_b
  end type octree_type
  
  interface octree_intersection_finder
    module procedure octree_intersection_finder_intersections, &
      & octree_intersection_finder_csr_sparsity
  end interface octree_intersection_finder
  
  interface allocate
    module procedure allocate_octree, allocate_node
  end interface allocate
  
  interface query
    module procedure query_node_internal, query_octree_allocatable, &
      & query_octree_pointer
  end interface query
  
  interface deallocate
    module procedure deallocate_octree, deallocate_node
  end interface deallocate
    
contains

  subroutine octree_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab, max_size)
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
    integer, optional, intent(in) :: max_size
    
    integer :: ele_a, nelements_a
    type(octree_type) :: octree
    
    if(size(positions_a, 1) /= TREE_DIM) then
      libsupermesh_abort("Invalid dimension")
    end if
    nelements_a = size(enlist_a, 2)
    if(nelements_a == 0 .or. size(enlist_b, 2) == 0) then
      do ele_a = 1, nelements_a
        allocate(map_ab(ele_a)%v(0))
        map_ab(ele_a)%n = 0
      end do
      return
    end if
    call allocate(octree, positions_b, enlist_b, max_size = max_size)
    
    do ele_a = 1, nelements_a
      call query(octree, positions_a(:, enlist_a(:, ele_a)), map_ab(ele_a)%v)
      map_ab(ele_a)%n = size(map_ab(ele_a)%v)
    end do
    
    call deallocate(octree)
    
  end subroutine octree_intersection_finder_intersections
  
  subroutine octree_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr, max_size)
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
    integer, optional, intent(in) :: max_size
    
    type(intersections), dimension(:), allocatable :: map_ab
    
    allocate(map_ab(size(enlist_a, 2)))
    call octree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab, max_size = max_size)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)
    
  end subroutine octree_intersection_finder_csr_sparsity

  subroutine allocate_node(root_node, positions, enlist, max_size)
    type(octree_node), intent(out) :: root_node
    ! TREE_DIM x nnodes
    real(kind = real_kind), dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    integer, optional, intent(in) :: max_size    
    
    integer, parameter :: default_max_size = 256 
    
    integer :: ele, nelements, node, nnodes
    
    if(size(positions, 1) /= TREE_DIM) then
      libsupermesh_abort("Invalid dimension")
    end if
    nnodes = size(positions, 2)
    if(nnodes == 0) then
      libsupermesh_abort("No nodes")
    end if
    if(size(enlist, 1) == 0) then
      libsupermesh_abort("No element local nodes")
    end if
    nelements = size(enlist, 2)
    
    if(present(max_size)) then
      ! A maximum number of stored elements has been provided. It is the
      ! caller's responsibility to ensure that this is large enough.
      root_node%max_n = max_size
    else
      root_node%max_n = max(max_nelist_degree(nnodes, enlist), default_max_size)
    end if
    allocate(root_node%elements(root_node%max_n))
    
    ! Compute the root node bounding box
    root_node%bbox(1, :) = positions(:, 1)
    root_node%bbox(2, :) = positions(:, 1)
    do node = 2, nnodes
      root_node%bbox(1, 1) = min(root_node%bbox(1, 1), positions(1, node))
      root_node%bbox(1, 2) = min(root_node%bbox(1, 2), positions(2, node))
      root_node%bbox(1, 3) = min(root_node%bbox(1, 3), positions(3, node))
      root_node%bbox(2, 1) = max(root_node%bbox(2, 1), positions(1, node))
      root_node%bbox(2, 2) = max(root_node%bbox(2, 2), positions(2, node))
      root_node%bbox(2, 3) = max(root_node%bbox(2, 3), positions(3, node))
    end do    
    
    ! Add all elements
    do ele = 1, nelements
      call add_element(root_node, ele, bbox(positions(:, enlist(:, ele))))
    end do
    
  end subroutine allocate_node
  
  subroutine allocate_octree(octree, positions, enlist, max_size)
    type(octree_type), intent(out) :: octree
    ! TREE_DIM x nnodes
    real(kind = real_kind), dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    integer, optional, intent(in) :: max_size  
    
    integer :: nelements
    
    nelements = size(enlist, 2)
    
    allocate(octree%octree)
    call allocate(octree%octree, positions, enlist, max_size = max_size)
    allocate(octree%eles_b(nelements), octree%neles_b, octree%seen_ele_b(nelements))
    octree%neles_b = 0
    octree%seen_ele_b = .false.
  
  end subroutine allocate_octree
  
  pure elemental subroutine deallocate_octree(octree)
    type(octree_type), intent(inout) :: octree
    
    deallocate(octree%eles_b, octree%neles_b, octree%seen_ele_b)
    call deallocate(octree%octree)
    deallocate(octree%octree)
  
  end subroutine deallocate_octree
  
  pure recursive subroutine deallocate_node(node)
    type(octree_node), intent(inout) :: node
    
    integer :: i
    
    if(node%n <= node%max_n) then
      deallocate(node%elements)
    else
      do i = 1, TREE_NCHILDREN
        call deallocate(node%children(i))
      end do
      deallocate(node%children)
      nullify(node%children)
    end if
    node%n = 0
    
  end subroutine deallocate_node
  
  pure recursive subroutine add_element(node, ele, bbox)
    type(octree_node), intent(inout) :: node
    integer, intent(in) :: ele
    real(kind = real_kind), dimension(2, TREE_DIM), intent(in) :: bbox
  
    integer :: i, j
    real(kind = real_kind), dimension(TREE_DIM) :: mid_point
  
    if(.not. bboxes_intersect(node%bbox, bbox)) return
    
    if(node%n <= node%max_n - 1) then
      ! Add new element to the leaf node
      node%n = node%n + 1
      node%elements(node%n)%bbox = bbox
      node%elements(node%n)%ele = ele
    else if(node%n == node%max_n) then
      ! Current leaf node is full. Split and add elements to the new child leaf
      ! nodes.
    
      ! Create new child leaf nodes
      allocate(node%children(TREE_NCHILDREN))
      mid_point = 0.5_real_kind * (node%bbox(1, :) + node%bbox(2, :))
      node%children(1)%bbox(1, :) = (/node%bbox(1, 1), node%bbox(1, 2), node%bbox(1, 3)/)
      node%children(1)%bbox(2, :) = (/mid_point(1),    mid_point(2),    mid_point(3)/)      
      node%children(2)%bbox(1, :) = (/node%bbox(1, 1), mid_point(2),    node%bbox(1, 3)/)
      node%children(2)%bbox(2, :) = (/mid_point(1),    node%bbox(2, 2), mid_point(3)/)
      node%children(3)%bbox(1, :) = (/mid_point(1),    mid_point(2),    node%bbox(1, 3)/)
      node%children(3)%bbox(2, :) = (/node%bbox(2, 1), node%bbox(2, 2), mid_point(3)/)
      node%children(4)%bbox(1, :) = (/mid_point(1),    node%bbox(1, 2), node%bbox(1, 3)/)
      node%children(4)%bbox(2, :) = (/node%bbox(2, 1), mid_point(2),    mid_point(3)/)
      node%children(5)%bbox(1, :) = (/node%bbox(1, 1), node%bbox(1, 2), mid_point(3)/)
      node%children(5)%bbox(2, :) = (/mid_point(1),    mid_point(2),    node%bbox(2, 3)/)      
      node%children(6)%bbox(1, :) = (/node%bbox(1, 1), mid_point(2),    mid_point(3)/)
      node%children(6)%bbox(2, :) = (/mid_point(1),    node%bbox(2, 2), node%bbox(2, 3)/)
      node%children(7)%bbox(1, :) = (/mid_point(1),    mid_point(2),    mid_point(3)/)
      node%children(7)%bbox(2, :) = (/node%bbox(2, 1), node%bbox(2, 2), node%bbox(2, 3)/)
      node%children(8)%bbox(1, :) = (/mid_point(1),    node%bbox(1, 2), mid_point(3)/)
      node%children(8)%bbox(2, :) = (/node%bbox(2, 1), mid_point(2),    node%bbox(2, 3)/)
      do i = 1, TREE_NCHILDREN
        node%children(i)%max_n = node%max_n
        allocate(node%children(i)%elements(node%max_n))
      end do
      
      ! Add elements stored in the parent to the new child leaf nodes
      do i = 1, node%max_n
        do j = 1, TREE_NCHILDREN
          call add_element(node%children(j), node%elements(i)%ele, node%elements(i)%bbox)
        end do
      end do
      ! Add new element to the new child leaf nodes
      do j = 1, TREE_NCHILDREN
        call add_element(node%children(j), ele, bbox)
      end do
      
      ! Mark the parent as not a leaf node
      node%n = huge(node%n)
      deallocate(node%elements)
      nullify(node%elements)
    else
      ! Not a leaf node. Add to child nodes.
      do j = 1, TREE_NCHILDREN
        call add_element(node%children(j), ele, bbox)
      end do
    end if
  
  end subroutine add_element
  
  pure recursive subroutine query_node_internal(node, bbox_a, eles_b, neles_b, seen_ele_b)
    type(octree_node), intent(in) :: node
    real(kind = real_kind), dimension(2, TREE_DIM), intent(in) :: bbox_a
    integer, dimension(:), intent(inout) :: eles_b
    integer, intent(inout) :: neles_b
    logical, dimension(:), intent(inout) :: seen_ele_b
    
    integer :: ele_b, i
    
    if(node%n == 0) return
    if(.not. bboxes_intersect(node%bbox, bbox_a)) return
    
    if(node%n <= node%max_n) then
      ! A leaf node. Return stored elements whose bounding boxes intersect with
      ! the query bounding box.
      do i = 1, node%n
        ele_b = node%elements(i)%ele
        if(.not. seen_ele_b(ele_b)) then
          if(bboxes_intersect(node%elements(i)%bbox, bbox_a)) then
            neles_b = neles_b + 1
            eles_b(neles_b) = ele_b
            seen_ele_b(ele_b) = .true.
          end if
        end if
      end do
    else
      ! Not a leaf node. Query child nodes.
      do i = 1, TREE_NCHILDREN
        call query(node%children(i), bbox_a, eles_b, neles_b, seen_ele_b)
      end do
    end if
    
  end subroutine query_node_internal
  
  pure subroutine query_octree_allocatable(octree, element_a, eles_b)
    type(octree_type), intent(inout) :: octree
    ! TREE_DIM x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    integer, dimension(:), allocatable, intent(out) :: eles_b
    
    octree%seen_ele_b(octree%eles_b(:octree%neles_b)) = .false.
    octree%neles_b = 0
    call query(octree%octree, bbox(element_a), octree%eles_b, octree%neles_b, octree%seen_ele_b)
    allocate(eles_b(octree%neles_b))
    eles_b = octree%eles_b(:octree%neles_b)
  
  end subroutine query_octree_allocatable
  
  pure subroutine query_octree_pointer(octree, element_a, eles_b)
    type(octree_type), intent(inout) :: octree
    ! TREE_DIM x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    integer, dimension(:), pointer, intent(out) :: eles_b
    
    octree%seen_ele_b(octree%eles_b(:octree%neles_b)) = .false.
    octree%neles_b = 0
    call query(octree%octree, bbox(element_a), octree%eles_b, octree%neles_b, octree%seen_ele_b)
    allocate(eles_b(octree%neles_b))
    eles_b = octree%eles_b(:octree%neles_b)
  
  end subroutine query_octree_pointer

  pure function bbox(coords)
    ! TREE_DIM x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: coords

    real(kind = real_kind), dimension(2, TREE_DIM) :: bbox

    integer :: i

    bbox(1, :) = coords(:, 1)
    bbox(2, :) = coords(:, 1)
    do i = 2, size(coords, 2)
      bbox(1, 1) = min(bbox(1, 1), coords(1, i))
      bbox(1, 2) = min(bbox(1, 2), coords(2, i))
      bbox(1, 3) = min(bbox(1, 3), coords(3, i))
      bbox(2, 1) = max(bbox(2, 1), coords(1, i))
      bbox(2, 2) = max(bbox(2, 2), coords(2, i))
      bbox(2, 3) = max(bbox(2, 3), coords(3, i))
    end do

  end function bbox

  pure function bboxes_intersect(bbox_1, bbox_2) result(intersect)
    real(kind = real_kind), dimension(2, TREE_DIM), intent(in) :: bbox_1
    real(kind = real_kind), dimension(2, TREE_DIM), intent(in) :: bbox_2

    logical :: intersect

    intersect = bbox_2(2, 1) > bbox_1(1, 1) .and. bbox_2(1, 1) < bbox_1(2, 1) &
        & .and. bbox_2(2, 2) > bbox_1(1, 2) .and. bbox_2(1, 2) < bbox_1(2, 2) &
        & .and. bbox_2(2, 3) > bbox_1(1, 3) .and. bbox_2(1, 3) < bbox_1(2, 3)

  end function bboxes_intersect

end module libsupermesh_octree_intersection_finder
