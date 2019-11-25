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
! femtools/tests/test_intersection_finder_completeness_3d.F90 in Fluidity git revision
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

subroutine test_intersection_finder_completeness_3d() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder, advancing_front_intersection_finder, &
    & octree_intersection_finder, tree_intersection_finder, &
    & rtree_intersection_finder, brute_force_intersection_finder
  use libsupermesh_precision, only : real_kind
  use libsupermesh_supermesh, only : tetrahedron_volume
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_tet_intersection, only : intersect_tets, tet_buf_size
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)

  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_tets_c
  real(kind = real_kind) :: volume_b, volume_c
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  logical :: fail
  real(kind = real_kind), dimension(3, 4, tet_buf_size) ::  tets_c
  real(kind = real_kind), dimension(:, :), allocatable :: positions_a, &
    & positions_b
  type(intersections), dimension(:), allocatable :: map_ba

  integer, parameter :: dim = 3

  call read_node("data/cube.1.node", dim, positions_a)
  call read_ele("data/cube.1.ele", dim, enlist_a)
  call read_node("data/cube.2.node", dim, positions_b)
  call read_ele("data/cube.2.ele", dim, enlist_b)
  allocate(map_ba(size(enlist_b, 2)))
  
  call intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("intersection_finder")
  call deallocate(map_ba)
  
  call advancing_front_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("advancing_front_intersection_finder")
  call deallocate(map_ba)
  
  call octree_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("octree_intersection_finder")
  call deallocate(map_ba)
  
  call rtree_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("rtree_intersection_finder")
  call deallocate(map_ba)
  
  call tree_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("tree_intersection_finder")
  call deallocate(map_ba)
  
  call brute_force_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("brute_force_intersection_finder")
  call deallocate(map_ba)

  deallocate(map_ba, positions_a, enlist_a, positions_b, enlist_b)

contains

  subroutine check_map_ba(name)
    character(len = *), intent(in) :: name

    fail = .false.
    do ele_b = 1, size(enlist_b, 2)
      volume_b = tetrahedron_volume(positions_b(:, enlist_b(:, ele_b)))

      volume_c = 0.0_real_kind
      do i = 1, map_ba(ele_b)%n
        ele_a = map_ba(ele_b)%v(i)
        call intersect_tets(positions_a(:, enlist_a(:, ele_a)), positions_b(:, enlist_b(:, ele_b)), tets_c, n_tets_c)
        do ele_c = 1, n_tets_c
          volume_c = volume_c + tetrahedron_volume(tets_c(:, :, ele_c))
        end do
      end do

      fail = (volume_b .fne. volume_c)
      if(fail) exit
    end do
    call report_test("[" // trim(name) // ": completeness]", fail, .false., "Need to have the same volume")
  
  end subroutine check_map_ba

end subroutine test_intersection_finder_completeness_3d
