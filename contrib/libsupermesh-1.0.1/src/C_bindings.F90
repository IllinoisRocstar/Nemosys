#include "libsupermesh_debug.h"

module libsupermesh_c_interface

  use iso_c_binding, only : c_double, c_long, c_int, c_float
  use libsupermesh_graphs, only : eelist_type, mesh_nelist, mesh_eelist, &
    & deallocate
  use libsupermesh_integer_set, only : integer_set, allocate, deallocate, &
    & insert, key_count, fetch
  use libsupermesh_intersection_finder, only : tree_intersection_finder, sort_intersection_finder
  use libsupermesh_tri_intersection, only : tri_buf_size, intersect_tris, triangle_area
  use libsupermesh_tet_intersection, only : tet_buf_size, intersect_tets, tetrahedron_volume
  use libsupermesh_interval_intersection, only : intersect_intervals, interval_size
  use libsupermesh_supermesh, only : intersect_tri_quad, intersect_tet_hex
  use libsupermesh_precision, only : real_kind

  implicit none

  private

  ! Compressed Sparse Row (CSR) sparsity patterns as described in:
  !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
  integer, dimension(:), allocatable, save :: nelist_indices
  ! nnodes + 1
  integer, dimension(:), allocatable, save :: nelist_indptr
  integer(kind = c_long), dimension(:), allocatable, save :: nnlist_indices
  ! nnodes + 1
  integer(kind = c_long), dimension(:), allocatable, save :: nnlist_indptr
  type(eelist_type), save :: eelist
  integer, dimension(:), allocatable, save :: map_ab_indices
  ! nelements_a + 1
  integer, dimension(:), allocatable, save :: map_ab_indptr

#ifdef LIBSUPERMESH_DOUBLE_PRECISION
   integer, parameter :: c_kind = c_double
#else
   integer, parameter :: c_kind = c_float
#endif

contains

  subroutine intersect_tet_hex_c( &
    & tet_a, hex_b, tets_c, n_tets_c) &
    & bind(c, name = "libsupermesh_intersect_tet_hex")
    real(kind = c_kind), dimension(3, 4), intent(in) :: tet_a
    real(kind = c_kind), dimension(3, 8), intent(in) :: hex_b
    real(kind = c_kind), dimension(3, 4, 729), intent(out) :: tets_c
    integer(kind = c_int), intent(out) :: n_tets_c

    call intersect_tet_hex(tet_a, hex_b, tets_c, n_tets_c) 

  end subroutine intersect_tet_hex_c

  subroutine intersect_tri_quad_c( &
    & tri_a, quad_b, tris_c, n_tris_c) &
    & bind(c, name = "libsupermesh_intersect_tri_quad")
    real(kind = c_kind), dimension(2, 3), intent(in) :: tri_a
    real(kind = c_kind), dimension(2, 4), intent(in) :: quad_b
    real(kind = c_kind), dimension(2, 3, 46), intent(out) :: tris_c
    integer(kind = c_int), intent(out) :: n_tris_c

    call intersect_tri_quad(tri_a, quad_b, tris_c, n_tris_c)

  end subroutine intersect_tri_quad_c

  subroutine sort_intersection_finder_set_input_c( &
    & nnodes_a, dim_a, nelements_a, loc_a, nnodes_b, dim_b, nelements_b, loc_b, &
    & positions_a, enlist_a, positions_b, enlist_b) &
    & bind(c, name = "libsupermesh_sort_intersection_finder_set_input")
    integer(kind = c_long), intent(in) :: nnodes_a
    integer(kind = c_int),  intent(in) :: dim_a
    integer(kind = c_long), intent(in) :: nelements_a
    integer(kind = c_int),  intent(in) :: loc_a
    integer(kind = c_long), intent(in) :: nnodes_b
    integer(kind = c_int),  intent(in) :: dim_b
    integer(kind = c_long), intent(in) :: nelements_b
    integer(kind = c_int),  intent(in) :: loc_b
    real(kind = c_kind), dimension(dim_a, nnodes_a), intent(in) :: positions_a
    integer(kind = c_long), dimension(loc_a, nelements_a), intent(in) :: enlist_a
    real(kind = c_kind), dimension(dim_b, nnodes_b), intent(in) :: positions_b
    integer(kind = c_long), dimension(loc_b, nelements_b), intent(in) :: enlist_b

    allocate(map_ab_indptr(nelements_a + 1))
    call sort_intersection_finder(real(positions_a, kind = real_kind), int(enlist_a + 1), real(positions_b, kind = real_kind), int(enlist_b + 1), map_ab_indices, map_ab_indptr)

  end subroutine sort_intersection_finder_set_input_c

  subroutine sort_intersection_finder_query_output_c( &
    & nindices) &
    bind(c, name = "libsupermesh_sort_intersection_finder_query_output")
    integer(kind = c_long), intent(out) :: nindices

    nindices = size(map_ab_indices)

  end subroutine sort_intersection_finder_query_output_c

  subroutine sort_intersection_finder_get_output_c( &
    & nelements, nindices, &
    & indices, ind_ptr) &
    & bind(c, name = "libsupermesh_sort_intersection_finder_get_output")
    integer(kind = c_long), intent(in) :: nelements
    integer(kind = c_long), intent(in) :: nindices
    integer(kind = c_long), dimension(nindices), intent(out) :: indices
    integer(kind = c_long), dimension(nelements + 1), intent(out) :: ind_ptr

    indices = map_ab_indices - 1
    ind_ptr = map_ab_indptr - 1

    deallocate(map_ab_indices, map_ab_indptr)

  end subroutine sort_intersection_finder_get_output_c

  subroutine intersect_intervals_c( &
    & interval_a, interval_b, intervals_c, n_intervals_c) &
    & bind(c, name = "libsupermesh_intersect_intervals")
    real(kind = c_kind), dimension(2), intent(in) :: interval_a
    real(kind = c_kind), dimension(2), intent(in) :: interval_b
    real(kind = c_kind), dimension(2), intent(out) :: intervals_c
    integer(kind = c_int), intent(out) :: n_intervals_c

    call intersect_intervals(interval_a, interval_b, intervals_c, n_intervals_c)

  end subroutine intersect_intervals_c

  subroutine interval_size_c(interval, size) &
    & bind(c, name = "libsupermesh_interval_size")
    real(kind = c_kind), dimension(2), intent(in) :: interval
    real(kind = c_kind), intent(out) :: size 

    size = interval_size(interval)

  end subroutine interval_size_c

  subroutine tetrahedron_volume_c(tet, volume) &
    & bind(c, name = "libsupermesh_tetrahedron_volume")
    real(kind = c_kind), dimension(3, 4), intent(in) :: tet
    real(kind = c_kind), intent(out) :: volume 

    volume = tetrahedron_volume(tet)

  end subroutine tetrahedron_volume_c

  subroutine triangle_area_c(tri, area) &
    & bind(c, name = "libsupermesh_triangle_area")
    real(kind = c_kind), dimension(2, 3), intent(in) :: tri
    real(kind = c_kind), intent(out) :: area 

    area = triangle_area(tri)

  end subroutine triangle_area_c

  subroutine intersect_tets_real_c( &
    & tet_a, tet_b, tets_c, n_tets_c) &
    & bind(c, name = "libsupermesh_intersect_tets_real")
    real(kind = c_kind), dimension(3, 4), intent(in) :: tet_a
    real(kind = c_kind), dimension(3, 4), intent(in) :: tet_b
    real(kind = c_kind), dimension(3, 4, tet_buf_size), intent(out) :: tets_c
    integer(kind = c_int), intent(out) :: n_tets_c

    call intersect_tets(tet_a, tet_b, tets_c, n_tets_c)

  end subroutine intersect_tets_real_c

  subroutine intersect_tris_real_c( &
    & tri_a, tri_b, tris_c, n_tris_c) &
    & bind(c, name = "libsupermesh_intersect_tris_real")
    real(kind = c_kind), dimension(2, 3), intent(in) :: tri_a
    real(kind = c_kind), dimension(2, 3), intent(in) :: tri_b
    real(kind = c_kind), dimension(2, 3, tri_buf_size), intent(out) :: tris_c
    integer(kind = c_int), intent(out) :: n_tris_c

    call intersect_tris(tri_a, tri_b, tris_c, n_tris_c)

  end subroutine intersect_tris_real_c

  subroutine mesh_nelist_set_input_c( &
    & nnodes, nelements, loc, &
    & enlist) &
    & bind(c, name = "libsupermesh_mesh_nelist_set_input")
    integer(kind = c_long), intent(in) :: nnodes
    integer(kind = c_long), intent(in) :: nelements
    integer(kind = c_long), intent(in) :: loc
    integer(kind = c_long), dimension(loc, nelements), intent(in) :: enlist

    allocate(nelist_indptr(nnodes + 1))
    call mesh_nelist(int(nnodes), int(enlist + 1), nelist_indices, nelist_indptr)

  end subroutine mesh_nelist_set_input_c

  subroutine mesh_nelist_query_output_c( &
    & nindices) &
    & bind(c, name = "libsupermesh_mesh_nelist_query_output")
    integer(kind = c_long), intent(out) :: nindices

    nindices = size(nelist_indices)

  end subroutine mesh_nelist_query_output_c

  subroutine mesh_nelist_get_output_c( &
    & nnodes, nindices, &
    & indices, ind_ptr) &
    & bind(c, name = "libsupermesh_mesh_nelist_get_output")
    integer(kind = c_long), intent(in) :: nnodes
    integer(kind = c_long), intent(in) :: nindices
    integer(kind = c_long), dimension(nindices), intent(out) :: indices
    integer(kind = c_long), dimension(nnodes + 1), intent(out) :: ind_ptr

    indices = nelist_indices - 1
    ind_ptr = nelist_indptr - 1

    deallocate(nelist_indices, nelist_indptr)

  end subroutine mesh_nelist_get_output_c

  subroutine mesh_nnlist(nnodes, enlist, nnlist_indices, nnlist_indptr)
    integer(kind = c_long), intent(in) :: nnodes
    ! loc x nelements
    integer(kind = c_long), dimension(:, :), intent(in) :: enlist
    integer(kind = c_long), dimension(:), allocatable, intent(out) :: nnlist_indices
    ! nnodes + 1
    integer(kind = c_long), dimension(:), intent(out) :: nnlist_indptr

    integer(kind = c_long) :: ele, loc, i, i_0, i_1, j, nelements, node
    type(integer_set), dimension(:), allocatable :: nnlist

    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    allocate(nnlist(nnodes))
    do node = 1, nnodes
      call allocate(nnlist(node))
    end do
    do ele = 1, nelements
      do i = 1, loc
        do j = i + 1, loc
          call insert(nnlist(enlist(i, ele)), int(enlist(j, ele)))
          call insert(nnlist(enlist(j, ele)), int(enlist(i, ele)))
        end do
      end do
    end do
    nnlist_indptr(1) = 1
    do node = 1, nnodes
      nnlist_indptr(node + 1) = nnlist_indptr(node) + key_count(nnlist(node))
    end do
    allocate(nnlist_indices(nnlist_indptr(nnodes + 1) - 1))
    do node = 1, nnodes
      i_0 = nnlist_indptr(node)
      i_1 = nnlist_indptr(node + 1) - 1
      do i = i_0, i_1
        nnlist_indices(i) = fetch(nnlist(node), int(i - i_0 + 1))
      end do
      call deallocate(nnlist(node))
    end do
    deallocate(nnlist)    

  end subroutine mesh_nnlist

  subroutine mesh_nnlist_set_input_c( &
    & nnodes, nelements, loc, &
    & enlist) &
    & bind(c, name = "libsupermesh_mesh_nnlist_set_input")
    integer(kind = c_long), intent(in) :: nnodes
    integer(kind = c_long), intent(in) :: nelements
    integer(kind = c_long), intent(in) :: loc
    integer(kind = c_long), dimension(loc, nelements), intent(in) :: enlist

    allocate(nnlist_indptr(nnodes + 1))
    call mesh_nnlist(nnodes, enlist + 1, nnlist_indices, nnlist_indptr)

  end subroutine mesh_nnlist_set_input_c

  subroutine mesh_nnlist_query_output_c( &
    & nindices) &
    & bind(c, name = "libsupermesh_mesh_nnlist_query_output")
    integer(kind = c_long), intent(out) :: nindices

    nindices = size(nnlist_indices)

  end subroutine mesh_nnlist_query_output_c

  subroutine mesh_nnlist_get_output_c( &
    & nnodes, nindices, &
    & indices, ind_ptr) &
    & bind(c, name = "libsupermesh_mesh_nnlist_get_output")
    integer(kind = c_long), intent(in) :: nnodes
    integer(kind = c_long), intent(in) :: nindices
    integer(kind = c_long), dimension(nindices), intent(out) :: indices
    integer(kind = c_long), dimension(nnodes + 1), intent(out) :: ind_ptr

    indices = nnlist_indices - 1
    ind_ptr = nnlist_indptr - 1

    deallocate(nnlist_indices, nnlist_indptr)

  end subroutine mesh_nnlist_get_output_c

  subroutine mesh_eelist_set_input_c( &
    & nnodes, nelements, loc, sloc, &
    & enlist) &
    & bind(c, name = "libsupermesh_mesh_eelist_set_input")
    integer(kind = c_long), intent(in) :: nnodes
    integer(kind = c_long), intent(in) :: nelements
    integer(kind = c_long), intent(in) :: loc
    integer(kind = c_long), intent(in) :: sloc
    integer(kind = c_long), dimension(loc, nelements), intent(in) :: enlist

    call mesh_eelist(int(nnodes), int(enlist + 1), int(sloc), eelist)

  end subroutine mesh_eelist_set_input_c

  subroutine mesh_eelist_query_output_c( &
    & nindices) &
    & bind(c, name = "libsupermesh_mesh_eelist_query_output")
    integer(kind = c_long), intent(out) :: nindices

    integer(kind = c_long) :: ele

    nindices = 0
    do ele = 1, size(eelist%n)
      nindices = nindices + eelist%n(ele)
    end do

  end subroutine mesh_eelist_query_output_c

  subroutine mesh_eelist_get_output_c( &
    & nelements, nindices, &
    & indices, ind_ptr) &
    & bind(c, name = "libsupermesh_mesh_eelist_get_output")
    integer(kind = c_long), intent(in) :: nelements
    integer(kind = c_long), intent(in) :: nindices
    integer(kind = c_long), dimension(nindices), intent(out) :: indices
    integer(kind = c_long), dimension(nelements + 1), intent(out) :: ind_ptr

    integer(kind = c_long) :: ele

    ind_ptr(1) = 0
    do ele = 1, nelements
      ind_ptr(ele + 1) = ind_ptr(ele) + eelist%n(ele)
      indices(ind_ptr(ele) + 1:ind_ptr(ele + 1)) = eelist%v(:eelist%n(ele), ele) - 1
    end do

    call deallocate(eelist)

  end subroutine mesh_eelist_get_output_c

  subroutine tree_intersection_finder_set_input_c( &
    & nnodes_a, dim_a, nelements_a, loc_a, nnodes_b, dim_b, nelements_b, loc_b, &
    & positions_a, enlist_a, positions_b, enlist_b) &
    & bind(c, name = "libsupermesh_tree_intersection_finder_set_input")
    integer(kind = c_long), intent(in) :: nnodes_a
    integer(kind = c_int),  intent(in) :: dim_a
    integer(kind = c_long), intent(in) :: nelements_a
    integer(kind = c_int),  intent(in) :: loc_a
    integer(kind = c_long), intent(in) :: nnodes_b
    integer(kind = c_int),  intent(in) :: dim_b
    integer(kind = c_long), intent(in) :: nelements_b
    integer(kind = c_int),  intent(in) :: loc_b
    real(kind = c_kind), dimension(dim_a, nnodes_a), intent(in) :: positions_a
    integer(kind = c_long), dimension(loc_a, nelements_a), intent(in) :: enlist_a
    real(kind = c_kind), dimension(dim_b, nnodes_b), intent(in) :: positions_b
    integer(kind = c_long), dimension(loc_b, nelements_b), intent(in) :: enlist_b

    allocate(map_ab_indptr(nelements_a + 1))
    call tree_intersection_finder(real(positions_a, kind = real_kind), int(enlist_a + 1), real(positions_b, kind = real_kind), int(enlist_b + 1), map_ab_indices, map_ab_indptr)

  end subroutine tree_intersection_finder_set_input_c

  subroutine tree_intersection_finder_query_output_c( &
    & nindices) &
    bind(c, name = "libsupermesh_tree_intersection_finder_query_output")
    integer(kind = c_long), intent(out) :: nindices

    nindices = size(map_ab_indices)

  end subroutine tree_intersection_finder_query_output_c

  subroutine tree_intersection_finder_get_output_c( &
    & nelements, nindices, &
    & indices, ind_ptr) &
    & bind(c, name = "libsupermesh_tree_intersection_finder_get_output")
    integer(kind = c_long), intent(in) :: nelements
    integer(kind = c_long), intent(in) :: nindices
    integer(kind = c_long), dimension(nindices), intent(out) :: indices
    integer(kind = c_long), dimension(nelements + 1), intent(out) :: ind_ptr

    indices = map_ab_indices - 1
    ind_ptr = map_ab_indptr - 1

    deallocate(map_ab_indices, map_ab_indptr)

  end subroutine tree_intersection_finder_get_output_c

  subroutine transpose_sparsity_c( &
    & nnodes, nnodes_t, nindices, &
    & indices, ind_ptr, indices_t, ind_ptr_t) &
    & bind(c, name = "libsupermesh_transpose_sparsity")
    integer(kind = c_long), intent(in) :: nnodes
    integer(kind = c_long), intent(in) :: nnodes_t
    integer(kind = c_long), intent(in) :: nindices
    integer(kind = c_long), dimension(nindices), intent(in) :: indices
    integer(kind = c_long), dimension(nnodes + 1), intent(in) :: ind_ptr
    integer(kind = c_long), dimension(nindices), intent(out) :: indices_t
    integer(kind = c_long), dimension(nnodes_t + 1), intent(out) :: ind_ptr_t

    integer(kind = c_long) :: i, node1, node2
    integer(kind = c_long), dimension(:), allocatable :: indices_n_t

    allocate(indices_n_t(nnodes_t))
    indices_n_t = 0
    do node1 = 1, nnodes
      do i = ind_ptr(node1) + 1, ind_ptr(node1 + 1)
        node2 = indices(i) + 1
        indices_n_t(node2) = indices_n_t(node2) + 1
      end do
    end do

    ind_ptr_t(1) = 0
    do node2 = 1, nnodes_t
      ind_ptr_t(node2 + 1) = ind_ptr_t(node2) + indices_n_t(node2)
    end do

    indices_n_t = 0
    do node1 = 1, nnodes
      do i = ind_ptr(node1) + 1, ind_ptr(node1 + 1)
        node2 = indices(i) + 1
        indices_t(ind_ptr_t(node2) + 1 + indices_n_t(node2)) = node1 - 1
        indices_n_t(node2) = indices_n_t(node2) + 1
      end do
    end do    
    deallocate(indices_n_t)

  end subroutine transpose_sparsity_c

end module libsupermesh_c_interface
