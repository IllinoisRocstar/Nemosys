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

subroutine test_parallel_supermesh_2d() bind(c)

  use iso_c_binding, only : c_int8_t
  use iso_fortran_env, only : error_unit, output_unit, real64
  use mpi

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_halo_ownership, only : element_ownership
  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_precision, only : mpi_real_kind, real_format, real_kind
  use libsupermesh_read_halos, only : halo_type, deallocate, read_halo
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_supermesh, only : tri_type, intersect_tris, tri_buf_size, &
    & triangle_area
  use libsupermesh_unittest_tools, only : operator(.fne.), report_test

  implicit none

  integer :: ele_a, ele_b, ele_c, i, ierr, loc_b, n_tris_c, nelements_a, &
    & nelements_b, nprocs, rank, real_extent
  integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real(kind = real_kind) :: area_parallel, area_serial, integral_parallel, &
    & integral_serial
  real(kind = real_kind), dimension(:), allocatable, target :: data_b_data, &
    & vals_b
  real(kind = real_kind), dimension(:, :), allocatable :: positions_a, &
    & positions_b
  type(halo_type) :: halo
  type(intersections), dimension(:), allocatable :: map_ab
  type(tri_type) :: tri_a, tri_b
  type(tri_type), dimension(tri_buf_size) :: tris_c
  
  character(len = int(log10(real(huge(rank), kind = real64))) + 2) :: rank_chr
  character(len = int(log10(real(huge(nprocs), kind = real64))) + 2) :: nprocs_chr

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr);  assert(ierr == MPI_SUCCESS)
  write(rank_chr, "(i0)") rank
  write(nprocs_chr, "(i0)") nprocs
  rank_chr = adjustl(rank_chr)
  nprocs_chr = adjustl(nprocs_chr)
  call MPI_Type_extent(mpi_real_kind, real_extent, ierr);  assert(ierr == MPI_SUCCESS)

  if(rank == 0) then
    call read_node("data/square_0_01.node", dim = 2, positions = positions_a)
    call read_ele("data/square_0_01.ele", dim = 2, enlist = enlist_a)
    call read_node("data/triangle_0_01.node", dim = 2, positions = positions_b)
    call read_ele("data/triangle_0_01.ele", dim = 2, enlist = enlist_b)
    nelements_a = size(enlist_a, 2)
    loc_b = size(enlist_b, 1)
    nelements_b = size(enlist_b, 2)
    allocate(vals_b(nelements_b))
    do ele_b = 1, nelements_b
      vals_b(ele_b) = sum(positions_b(1, enlist_b(:, ele_b))) / real(loc_b, kind = real_kind)
    end do
    
    allocate(map_ab(nelements_a))
    call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

    area_serial = 0.0_real_kind
    integral_serial = 0.0_real_kind
    do ele_a = 1, nelements_a
      tri_a%v = positions_a(:, enlist_a(:, ele_a))

      do i = 1, map_ab(ele_a)%n
        ele_b = map_ab(ele_a)%v(i)
        tri_b%v = positions_b(:, enlist_b(:, ele_b))

        call intersect_tris(tri_a, tri_b, tris_c, n_tris_c)

        do ele_c = 1, n_tris_c
          area_serial = area_serial + triangle_area(tris_c(ele_c))
          integral_serial = integral_serial + vals_b(ele_b) * triangle_area(tris_c(ele_c))
        end do
      end do
    end do

    call deallocate(map_ab)
    deallocate(map_ab)

    deallocate(positions_a, enlist_a, positions_b, enlist_b, vals_b)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  call read_node("data/square_0_01_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".node", dim = 2, positions = positions_a)
  call read_ele("data/square_0_01_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_a)
  call read_halo("data/square_0_01_" // trim(nprocs_chr), halo, level = 2)
  nelements_a = size(enlist_a, 2)
  allocate(ele_owner_a(nelements_a))
  call element_ownership(size(positions_a, 2), enlist_a, halo, ele_owner_a)
  call deallocate(halo)

  call read_node("data/triangle_0_01_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".node", dim = 2, positions = positions_b)
  call read_ele("data/triangle_0_01_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_b)
  call read_halo("data/triangle_0_01_" // trim(nprocs_chr), halo, level = 2)
  loc_b = size(enlist_b, 1)
  nelements_b = size(enlist_b, 2)
  allocate(ele_owner_b(nelements_b))
  call element_ownership(size(positions_b, 2), enlist_b, halo, ele_owner_b)
  call deallocate(halo)

  allocate(vals_b(nelements_b))
  do ele_b = 1, nelements_b
    vals_b(ele_b) = sum(positions_b(1, enlist_b(:, ele_b))) / real(loc_b, kind = real_kind)
  end do
  
  area_parallel = 0.0_real_kind
  integral_parallel = 0.0_real_kind
  call parallel_supermesh(positions_a, enlist_a, ele_owner_a, &
                        & positions_b, enlist_b, ele_owner_b, &
                        & pack_data_b, unpack_data_b, intersection_calculation)
  call cleanup_data_b()
  call MPI_Allreduce(MPI_IN_PLACE, area_parallel, 1, mpi_real_kind, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, mpi_real_kind, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  deallocate(positions_a, enlist_a, ele_owner_a)
  deallocate(positions_b, enlist_b, ele_owner_b, vals_b)

  flush(output_unit)
  flush(error_unit)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  if(rank == 0) then
    write(output_unit, "(a,"//real_format//")") "Area, serial       =", area_serial
    write(output_unit, "(a,"//real_format//")") "Area, parallel     =", area_parallel
    write(output_unit, "(a,"//real_format//")") "Integral, serial   =", integral_serial
    write(output_unit, "(a,"//real_format//")") "Integral, parallel =", integral_parallel

    call report_test("[test_parallel_supermesh serial area]", &
      & area_serial .fne. 0.5_real_kind, .false., "Incorrect serial area")
    call report_test("[test_parallel_supermesh areas]", &
      & area_serial .fne. area_parallel, .false., &
      & "Serial and parallel areas differ")
    call report_test("[test_parallel_supermesh integrals]", &
      & integral_serial .fne. integral_parallel, .false., &
      & "Serial and parallel integrals differ")
  end if
  flush(output_unit)
  flush(error_unit)

contains

  subroutine pack_data_b(nodes_b, eles_b, data_b)
    integer, dimension(:), intent(in) :: nodes_b
    integer, dimension(:), intent(in) :: eles_b
    integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_b

    integer :: ierr, ndata_b, position
    real(kind = real_kind), dimension(:), allocatable :: ldata

    allocate(ldata(size(eles_b)))
    ldata = vals_b(eles_b)

    ndata_b = size(eles_b) * real_extent
    allocate(data_b(ndata_b))
    position = 0
    call MPI_Pack(ldata, size(eles_b), mpi_real_kind, data_b, ndata_b, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

    deallocate(ldata)

  end subroutine pack_data_b

  subroutine unpack_data_b(nnodes_b, nelements_b, data_b)
    integer, intent(in) :: nnodes_b
    integer, intent(in) :: nelements_b
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    integer :: ierr, position

    call cleanup_data_b()

    allocate(data_b_data(nelements_b))
    position = 0
    call MPI_Unpack(data_b, size(data_b), position, data_b_data, size(data_b_data), mpi_real_kind, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  end subroutine unpack_data_b

  subroutine cleanup_data_b()
    if(allocated(data_b_data)) deallocate(data_b_data)
  end subroutine cleanup_data_b

  subroutine intersection_calculation(element_a, element_b, elements_c, nodes_b, ele_a, ele_b, local)
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    ! dim x loc_b
    real(kind = real_kind), dimension(:, :), intent(in) :: element_b
    ! dim x loc_c x nelements_c
    real(kind = real_kind), dimension(:, :, :), intent(in) :: elements_c
    ! loc_b
    integer, dimension(:), intent(in) :: nodes_b
    integer, intent(in) :: ele_a
    integer, intent(in) :: ele_b
    logical, intent(in) :: local

    integer :: ele_c, nelements_c
    real(kind = real_kind), dimension(:), pointer :: data_b
    
    if(local) then
      data_b => vals_b
    else
      data_b => data_b_data
    end if

    nelements_c = size(elements_c, 3)
    do ele_c = 1, nelements_c
      area_parallel = area_parallel + triangle_area(elements_c(:, :, ele_c))
      integral_parallel = integral_parallel + data_b(ele_b) * triangle_area(elements_c(:, :, ele_c))
    end do

  end subroutine intersection_calculation

end subroutine test_parallel_supermesh_2d
