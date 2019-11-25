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

subroutine test_parallel_p2_inner_product_2d() bind(c)

  use iso_c_binding, only : c_int8_t
  use iso_fortran_env, only : error_unit, output_unit, real64
  use mpi

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_halo_ownership, only : element_ownership
  use libsupermesh_integer_map, only : integer_map, allocate, deallocate, &
    & has_key, fetch, insert
  use libsupermesh_integer_set, only : integer_set, allocate, deallocate, &
    & insert, key_count, fetch
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_precision, only : mpi_real_kind, real_format, real_kind
  use libsupermesh_read_halos, only : halo_type, deallocate, read_halo
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_supermesh, only : triangle_area
  use libsupermesh_unittest_tools, only : operator(.fne.), report_test

  implicit none
  
  ! Input Triangle mesh base names
  character(len = *), parameter :: basename_a = "data/triangle_0_01", &
                                 & basename_b = "data/square_0_01"
  real(kind = real_kind), parameter :: area_ref = 0.5_real_kind, &
    & integral_ref = 2.7083333333333272e-02_real_kind

  integer :: ierr, integer_extent, nprocs, rank, real_extent

  integer :: nelements_a, nelements_b, nnodes_p1_a, nnodes_p2_a, &
    & nnodes_p1_b, nnodes_p2_b
  integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
  integer, dimension(:, :), allocatable :: enlist_p1_a, enlist_p1_b, &
    & enlist_p2_a
  integer, dimension(:, :), allocatable, target :: enlist_p2_b
  real(kind = real_kind), dimension(:), allocatable :: interpolated, field_a
  real(kind = real_kind), dimension(:), allocatable, target :: field_b
  real(kind = real_kind), dimension(:, :), allocatable :: positions_a, positions_b
  type(halo_type) :: halo_a, halo_b

  integer :: data_nelements_b, data_nnodes_p2_b
  integer, dimension(:, :), allocatable, target :: data_enlist_p2_b
  real(kind = real_kind), dimension(:), allocatable, target :: data_field_b

  ! P2 mass matrix
  real(kind = real_kind), dimension(6, 6), parameter :: mass_p2 = &
    & reshape((/ 6.0_real_kind, -1.0_real_kind, -1.0_real_kind,  0.0_real_kind, -4.0_real_kind,  0.0_real_kind, &
              & -1.0_real_kind,  6.0_real_kind, -1.0_real_kind,  0.0_real_kind,  0.0_real_kind, -4.0_real_kind, &
              & -1.0_real_kind, -1.0_real_kind,  6.0_real_kind, -4.0_real_kind,  0.0_real_kind,  0.0_real_kind, &
              &  0.0_real_kind,  0.0_real_kind, -4.0_real_kind, 32.0_real_kind, 16.0_real_kind, 16.0_real_kind, &
              & -4.0_real_kind,  0.0_real_kind,  0.0_real_kind, 16.0_real_kind, 32.0_real_kind, 16.0_real_kind, &
              &  0.0_real_kind, -4.0_real_kind,  0.0_real_kind, 16.0_real_kind, 16.0_real_kind, 32.0_real_kind/) &
    & / 360.0_real_kind, (/6, 6/))
  real(kind = real_kind) :: area_parallel, integral_parallel
  
  character(len = int(log10(real(huge(rank), kind = real64))) + 2) :: rank_chr
  character(len = int(log10(real(huge(nprocs), kind = real64))) + 2) :: nprocs_chr

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr);  assert(ierr == MPI_SUCCESS)
  write(rank_chr, "(i0)") rank
  write(nprocs_chr, "(i0)") nprocs
  rank_chr = adjustl(rank_chr)
  nprocs_chr = adjustl(nprocs_chr)
  call MPI_Type_extent(MPI_INTEGER, integer_extent, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Type_extent(mpi_real_kind, real_extent, ierr);  assert(ierr == MPI_SUCCESS)

  ! Read the donor mesh partition
  call read_node(trim(basename_a) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".node", dim = 2, positions = positions_a)
  call read_ele(trim(basename_a) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_p1_a)
  nnodes_p1_a = size(positions_a, 2)
  nelements_a = size(enlist_p1_a, 2)
  ! Read donor mesh halo data ...
  call read_halo(trim(basename_a) // "_" // trim(nprocs_chr), halo_a, level = 2)
  ! ... and determine the donor mesh element ownership
  allocate(ele_owner_a(nelements_a))
  call element_ownership(nnodes_p1_a, enlist_p1_a, halo_a, ele_owner_a)
  ! Generate the donor P2 element-node graph
  allocate(enlist_p2_a(6, nelements_a))
  call p2_connectivity(nnodes_p1_a, enlist_p1_a, nnodes_p2_a, enlist_p2_a)
  ! Construct a donor P2 field equal to: x y
  allocate(field_a(nnodes_p2_a), interpolated(nnodes_p2_a))
  call interpolate_p1_p2(enlist_p1_a, positions_a(1, :), enlist_p2_a, field_a)
  call interpolate_p1_p2(enlist_p1_a, positions_a(2, :), enlist_p2_a, interpolated)
  field_a = field_a * interpolated
  deallocate(interpolated)


  ! Read the target mesh partition
  call read_node(trim(basename_b) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".node", dim = 2, positions = positions_b)
  call read_ele(trim(basename_b) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_p1_b)
  nnodes_p1_b = size(positions_b, 2)
  nelements_b = size(enlist_p1_b, 2)
  ! Read target mesh halo data ...
  call read_halo(trim(basename_b) // "_" // trim(nprocs_chr), halo_b, level = 2)
  ! ... and determine the target mesh element ownership
  allocate(ele_owner_b(nelements_b))
  call element_ownership(nnodes_p1_b, enlist_p1_b, halo_b, ele_owner_b)
  ! Generate the donor P2 element-node graph
  allocate(enlist_p2_b(6, nelements_b))
  call p2_connectivity(nnodes_p1_b, enlist_p1_b, nnodes_p2_b, enlist_p2_b)
  ! Construct a target P2 field equal to: x^2
  allocate(field_b(nnodes_p2_b))
  call interpolate_p1_p2(enlist_p1_b, positions_b(1, :), enlist_p2_b, field_b)
  field_b = field_b * field_b

  ! Perform multi-mesh integration
  area_parallel = 0.0_real_kind
  integral_parallel = 0.0_real_kind
  call parallel_supermesh(positions_a, enlist_p1_a, ele_owner_a, &
                        & positions_b, enlist_p1_b, ele_owner_b, &
                        & pack_data_b, unpack_data_b, intersection_calculation, &
                        & comm = MPI_COMM_WORLD)
  ! Deallocate any remaining unpacked communicated data
  call cleanup_unpack_data_b()

  ! Sum all process contributions to the multi-mesh integrals
  call MPI_Allreduce(MPI_IN_PLACE, area_parallel, 1, mpi_real_kind, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, mpi_real_kind, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  flush(output_unit)
  flush(error_unit)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  if(rank == 0) then
    ! Display the multi-mesh integrals on rank 0
    write(output_unit, "(a,"//real_format//",a,"//real_format//",a)") &
      & "Area     = ", area_parallel, " (error = ", abs(area_parallel - area_ref), ")"
    write(output_unit, "(a,"//real_format//",a,"//real_format//",a)") &
      & "Integral = ", integral_parallel, " (error = ", abs(integral_parallel - integral_ref), ")"
    
    ! Test the multi-mesh integrals against the reference values
    call report_test("[test_parallel_p1_inner_product_2d area]", area_parallel .fne. area_ref, .false., "Incorrect area")
    call report_test("[test_parallel_p1_inner_product_2d integral]", integral_parallel .fne. integral_ref, .false., "Incorrect integral")
  end if
  flush(output_unit)
  flush(error_unit)

  ! Cleanup
  deallocate(positions_a, enlist_p1_a, ele_owner_a, enlist_p2_a, field_a, &
           & positions_b, enlist_p1_b, ele_owner_b, enlist_p2_b, field_b)
  call deallocate(halo_a)
  call deallocate(halo_b)

contains

  ! Generate a P2 element-node graph
  subroutine p2_connectivity(nnodes_p1, enlist_p1, nnodes_p2, enlist_p2)
    ! Number of P1 nodes
    integer, intent(in) :: nnodes_p1
    ! P1 element-node graph
    ! Shape: 3 x nelements
    integer, dimension(:, :), intent(in) :: enlist_p1
    ! Number of P2 nodes
    integer, intent(out) :: nnodes_p2
    ! P2 element-node graph
    ! Shape: 6 x nelements
    integer, dimension(:, :), intent(out) :: enlist_p2
    
    integer :: ele, lnode, nelements, node_p1, node_p1_1, node_p1_2, node_p2
    type(integer_map) :: node_map_p1_p2
    type(integer_map), dimension(:), allocatable :: node_map_p2
    
    nelements = size(enlist_p1, 2)
    
    call allocate(node_map_p1_p2)
    allocate(node_map_p2(nnodes_p1))
    call allocate(node_map_p2)
    nnodes_p2 = 0
        
    do ele = 1, nelements    
      do lnode = 1, 6
        select case(lnode)
          case(1:3)
            node_p1 = enlist_p1(lnode, ele)
            if(has_key(node_map_p1_p2, node_p1)) then
              node_p2 = fetch(node_map_p1_p2, node_p1)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p1_p2, node_p1, node_p2)
            end if
          case(4)
            node_p1_1 = enlist_p1(1, ele)
            node_p1_2 = enlist_p1(2, ele)
            node_p1 = min(node_p1_1, node_p1_2)
            node_p1_2 = max(node_p1_1, node_p1_2)
            node_p1_1 = node_p1
            if(has_key(node_map_p2(node_p1_1), node_p1_2)) then
              node_p2 = fetch(node_map_p2(node_p1_1), node_p1_2)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p2(node_p1_1), node_p1_2, node_p2)
            end if
          case(5)
            node_p1_1 = enlist_p1(2, ele)
            node_p1_2 = enlist_p1(3, ele)
            node_p1 = min(node_p1_1, node_p1_2)
            node_p1_2 = max(node_p1_1, node_p1_2)
            node_p1_1 = node_p1
            if(has_key(node_map_p2(node_p1_1), node_p1_2)) then
              node_p2 = fetch(node_map_p2(node_p1_1), node_p1_2)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p2(node_p1_1), node_p1_2, node_p2)
            end if
          case(6)
            node_p1_1 = enlist_p1(1, ele)
            node_p1_2 = enlist_p1(3, ele)
            node_p1 = min(node_p1_1, node_p1_2)
            node_p1_2 = max(node_p1_1, node_p1_2)
            node_p1_1 = node_p1
            if(has_key(node_map_p2(node_p1_1), node_p1_2)) then
              node_p2 = fetch(node_map_p2(node_p1_1), node_p1_2)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p2(node_p1_1), node_p1_2, node_p2)
            end if
          case default
            libsupermesh_abort("Invalid local node")
        end select
        enlist_p2(lnode, ele) = node_p2
      end do      
    end do
    
    call deallocate(node_map_p1_p2)
    call deallocate(node_map_p2)
    deallocate(node_map_p2)
    
  end subroutine p2_connectivity
  
  ! Interpolate a P1 field onto a P2 function space
  subroutine interpolate_p1_p2(enlist_p1, field_p1, enlist_p2, field_p2)
    ! P1 element-node graph
    ! Shape: 3 x nelements
    integer, dimension(:, :), intent(in) :: enlist_p1
    ! P1 field
    ! Shape: nnodes_p1
    real(kind = real_kind), dimension(:), intent(in) :: field_p1
    ! P2 element-node graph
    ! Shape: 6 x nelements
    integer, dimension(:, :), intent(in) :: enlist_p2
    ! P2 field
    ! Shape: nnodes_p2
    real(kind = real_kind), dimension(:), intent(out) :: field_p2
    
    integer :: ele, lnode, nelements, node_p2, nnodes_p2
    logical, dimension(:), allocatable :: seen_node_p2
    
    nelements = size(enlist_p2, 2)
    nnodes_p2 = size(field_p2)
    
    allocate(seen_node_p2(nnodes_p2))
    seen_node_p2 = .false.
    
    do ele = 1, nelements
      do lnode = 1, 6
        node_p2 = enlist_p2(lnode, ele)
        if(seen_node_p2(node_p2)) cycle
        seen_node_p2(node_p2) = .true.
        select case(lnode)
          case(1:3)
            field_p2(node_p2) = field_p1(enlist_p1(lnode, ele))
          case(4)
            field_p2(node_p2) = 0.5_real_kind * (field_p1(enlist_p1(1, ele)) + field_p1(enlist_p1(2, ele)))
          case(5)
            field_p2(node_p2) = 0.5_real_kind * (field_p1(enlist_p1(2, ele)) + field_p1(enlist_p1(3, ele)))
          case(6)
            field_p2(node_p2) = 0.5_real_kind * (field_p1(enlist_p1(1, ele)) + field_p1(enlist_p1(3, ele)))
          case default
            libsupermesh_abort("Invalid local node")
        end select
      end do
    end do
    
    deallocate(seen_node_p2)
  
  end subroutine interpolate_p1_p2

  ! Given the provided mesh vertices and elements, pack data for communication
  subroutine pack_data_b(nodes_b, eles_b, data_b)
    ! Mesh vertices to be communicated
    integer, dimension(:), intent(in) :: nodes_b
    ! Mesh elements to be communicated
    integer, dimension(:), intent(in) :: eles_b
    ! Packed data for communication
    integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_b
 
    integer :: ele, i, lnode, ndata_b, node_b, nnodes_p2_b, position
    integer, dimension(:, :), allocatable :: data_enlist_p2_b
    type(integer_map) :: node_map
    type(integer_set) :: nodes_p2_b
    
    ! For which P2 nodes do we need to send data?
    call allocate(nodes_p2_b)
    do i = 1, size(eles_b)
      call insert(nodes_p2_b, enlist_p2_b(:, eles_b(i)))
    end do
    nnodes_p2_b = key_count(nodes_p2_b)
    
    ! Construct a map from communicated nodes to local nodes ...
    call allocate(node_map)
    do i = 1, nnodes_p2_b
      node_b = fetch(nodes_p2_b, i)
      call insert(node_map, node_b, i)
    end do
    ! ... and use this to construct the communicated P2 element-node graph
    allocate(data_enlist_p2_b(6, size(eles_b)))
    do i = 1, size(eles_b)
      ele = eles_b(i)
      do lnode = 1, 6
        node_b = enlist_p2_b(lnode, ele)
        data_enlist_p2_b(lnode, i) = fetch(node_map, node_b)
      end do
    end do
    
    ! Gather P2 field values for communication
    allocate(data_field_b(nnodes_p2_b))
    do i = 1, nnodes_p2_b
      node_b = fetch(nodes_p2_b, i)
      data_field_b(i) = field_b(node_b)
    end do
    
    ! Pack data for communication:
    !   1 integer                         -- number of P2 nodes
    !   (6 x number of elements) integers -- communicated P2 element-node graph
    !   (number of P2 nodes) reals        -- communicated P2 field values
    ndata_b = (1 + 6 * size(eles_b)) * integer_extent + nnodes_p2_b * real_extent
    allocate(data_b(ndata_b))
    position = 0
    call MPI_Pack(nnodes_p2_b, 1, MPI_INTEGER, data_b, ndata_b, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Pack(data_enlist_p2_b, 6 * size(eles_b), MPI_INTEGER, data_b, ndata_b, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Pack(data_field_b, nnodes_p2_b, mpi_real_kind, data_b, ndata_b, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    
    call deallocate(nodes_p2_b)
    call deallocate(node_map)
    deallocate(data_enlist_p2_b, data_field_b)
    
  end subroutine pack_data_b
  
  ! Unpack communicated data
  subroutine unpack_data_b(nnodes_b, nelements_b, data_b)
    ! Number of communicated mesh vertices
    integer, intent(in) :: nnodes_b
    ! Number of communicated elements
    integer, intent(in) :: nelements_b
    ! Packed communicated data
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b
    
    integer :: position
    
    ! Deallocate any previously unpacked communicated data
    call cleanup_unpack_data_b()
    
    position = 0
    ! Store the number of elements
    data_nelements_b = nelements_b
    ! Unpack the number of P2 nodes
    call MPI_Unpack(data_b, size(data_b), position, data_nnodes_p2_b, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    ! Unpack the P2 element-node graph
    allocate(data_enlist_p2_b(6, data_nelements_b))
    call MPI_Unpack(data_b, size(data_b), position, data_enlist_p2_b, 6 * data_nelements_b, MPI_INTEGER, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    ! Unpack the P2 field values
    allocate(data_field_b(data_nnodes_p2_b))
    call MPI_Unpack(data_b, size(data_b), position, data_field_b, data_nnodes_p2_b, mpi_real_kind, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    
  end subroutine unpack_data_b
  
  ! Deallocate any previously unpacked communicated data
  subroutine cleanup_unpack_data_b()
    if(allocated(data_enlist_p2_b)) deallocate(data_enlist_p2_b)
    if(allocated(data_field_b)) deallocate(data_field_b)
  end subroutine cleanup_unpack_data_b

  ! Evaluate P2 basis functions at a given point
  pure function basis_functions_p2(cell_coords, coord) result(fns)
    ! Triangle vertex coordinates
    ! Shape: dim x loc_p1
    real(kind = real_kind), dimension(2, 3), intent(in) :: cell_coords
    ! Coordinate at which to evaluate the basis functions
    ! Shape: dim
    real(kind = real_kind), dimension(2), intent(in) :: coord

    ! Basis function values
    ! Shape: loc_p2
    real(kind = real_kind), dimension(6) :: fns

    real(kind = real_kind) :: x, y
    real(kind = real_kind), dimension(2) :: e1, e2, lcoord
    real(kind = real_kind), dimension(2, 2) :: jac

    e1 = cell_coords(:, 2) - cell_coords(:, 1)
    e2 = cell_coords(:, 3) - cell_coords(:, 1)

    jac(1, 1) =  e2(2);  jac(1, 2) = -e2(1)
    jac(2, 1) = -e1(2);  jac(2, 2) =  e1(1)
    jac = jac / (jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1))

    lcoord = matmul(jac, coord - cell_coords(:, 1))
    x = lcoord(1)
    y = lcoord(2)
    
    fns(1) = (1.0_real_kind - x - y) * (1.0_real_kind - 2.0_real_kind * x - 2.0_real_kind * y)
    fns(2) = x * (2.0_real_kind * x - 1.0_real_kind)
    fns(3) = y * (2.0_real_kind * y - 1.0_real_kind)
    fns(4) = 4.0_real_kind * x * (1.0_real_kind - x - y)
    fns(5) = 4.0_real_kind * x * y
    fns(6) = 4.0_real_kind * y * (1.0_real_kind - x - y)

  end function basis_functions_p2

  ! Interpolate a P2 function at given point
  pure function interpolate_p2(cell_coords_d, cell_x_d, coord_s) result(x_s)
    ! Triangle vertex coordinates
    ! Shape: dim x loc_p1
    real(kind = real_kind), dimension(2, 3), intent(in) :: cell_coords_d
    ! P2 nodal values
    ! Shape: loc_p2
    real(kind = real_kind), dimension(6), intent(in) :: cell_x_d
    ! Coordinate at which to evaluate the P2 function
    ! Shape: dim
    real(kind = real_kind), dimension(2), intent(in) :: coord_s

    real(kind = real_kind) :: x_s

    x_s = dot_product(basis_functions_p2(cell_coords_d, coord_s), cell_x_d)

  end function interpolate_p2
  
  ! Perform calculations on the local supermesh
  subroutine intersection_calculation(element_a, element_b, elements_c, nodes_b, ele_a, ele_b, local)
    ! Target mesh element vertex coordinates
    ! Shape: dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    ! Donor mesh element vertex coordinates
    ! Shape: dim x loc_b
    real(kind = real_kind), dimension(:, :), intent(in) :: element_b
    ! Supermesh element vertex coordinates
    ! Shape: dim x loc_c x nelements_c
    real(kind = real_kind), dimension(:, :, :), intent(in) :: elements_c
    ! Donor mesh vertex indices
    ! Shape: loc_b
    integer, dimension(:), intent(in) :: nodes_b
    ! Target mesh element
    integer, intent(in) :: ele_a
    ! Donor mesh element
    integer, intent(in) :: ele_b
    ! Whether this is a local calculation or a calculation using communicated
    ! data
    logical, intent(in) :: local
    
    integer :: ele_c, lnode
    integer, dimension(:, :), pointer :: lenlist_p2_b
    real(kind = real_kind) :: area
    real(kind = real_kind), dimension(6) :: field_a_c, field_b_c
    real(kind = real_kind), dimension(:), pointer :: lfield_b
    
    if(local) then
      ! If this is a local calculation, use the local P2 element-node graph and
      ! field data
      lenlist_p2_b => enlist_p2_b
      lfield_b => field_b
    else
      ! Otherwise, use the unpacked communicated element-node graph and P2 field
      ! data
      lenlist_p2_b => data_enlist_p2_b
      lfield_b => data_field_b
    end if
  
    do ele_c = 1, size(elements_c, 3)
      ! Compute the supermesh triangle area
      area = triangle_area(elements_c(:, :, ele_c))
      ! Local contribution to the intersection area
      area_parallel = area_parallel + area
      ! Interpolate the donor and target P2 functions onto a P2 space on the
      ! supermesh element
      do lnode = 1, 6
        select case(lnode)
          case(1:3)
            field_b_c(lnode) = interpolate_p2(element_b, lfield_b(lenlist_p2_b(:, ele_b)), elements_c(:, lnode, ele_c))
            field_a_c(lnode) = interpolate_p2(element_a, field_a(enlist_p2_a(:, ele_a)), elements_c(:, lnode, ele_c))
          case(4)
            field_b_c(lnode) = interpolate_p2(element_b, lfield_b(lenlist_p2_b(:, ele_b)), 0.5_real_kind * (elements_c(:, 1, ele_c) + elements_c(:, 2, ele_c)))
            field_a_c(lnode) = interpolate_p2(element_a, field_a(enlist_p2_a(:, ele_a)), 0.5_real_kind * (elements_c(:, 1, ele_c) + elements_c(:, 2, ele_c)))
          case(5)
            field_b_c(lnode) = interpolate_p2(element_b, lfield_b(lenlist_p2_b(:, ele_b)), 0.5_real_kind * (elements_c(:, 2, ele_c) + elements_c(:, 3, ele_c)))
            field_a_c(lnode) = interpolate_p2(element_a, field_a(enlist_p2_a(:, ele_a)), 0.5_real_kind * (elements_c(:, 2, ele_c) + elements_c(:, 3, ele_c)))
          case(6)
            field_b_c(lnode) = interpolate_p2(element_b, lfield_b(lenlist_p2_b(:, ele_b)), 0.5_real_kind * (elements_c(:, 1, ele_c) + elements_c(:, 3, ele_c)))
            field_a_c(lnode) = interpolate_p2(element_a, field_a(enlist_p2_a(:, ele_a)), 0.5_real_kind * (elements_c(:, 1, ele_c) + elements_c(:, 3, ele_c)))
          case default
            libsupermesh_abort("Invalid local node")
        end select
      end do
      ! Local contribution to the multi-mesh inner product
      integral_parallel = integral_parallel + 2.0_real_kind * area * dot_product(field_b_c, matmul(mass_p2, field_a_c))
    end do
    
  end subroutine intersection_calculation
  
end subroutine test_parallel_p2_inner_product_2d
