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
! femtools/tests/test_tet_intersector.F90 in Fluidity git revision
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

subroutine test_tet_intersector() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_precision, only : real_kind, compensated_sum, initialise, &
    & add, sum
  use libsupermesh_supermesh, only : intersect_elements, &
    & intersect_simplices, tetrahedron_volume
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_tet_intersection, only : tet_type, tet_buf_size, &
    & intersect_tets, intersect_polys, get_planes
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)
  
  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_tets_c
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real(kind = real_kind), dimension(3, 4) :: tet_a_real
  real(kind = real_kind), dimension(3, 4, tet_buf_size) :: tets_c_real
  real(kind = real_kind), dimension(:, :), allocatable :: positions_a, &
    & positions_b 
  type(compensated_sum) :: volume_c
  type(intersections), dimension(:), allocatable :: map_ab
  type(tet_type) :: tet_a, tet_b
  type(tet_type), dimension(tet_buf_size) :: tets_c, work

  integer, parameter :: dim = 3
  
  call read_node("data/pyramid_0_05.node", dim, positions_a)
  call read_ele("data/pyramid_0_05.ele", dim, enlist_a)
  call read_node("data/cube_0_05.node", dim, positions_b)
  call read_ele("data/cube_0_05.ele", dim, enlist_b)
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  call initialise(volume_c)
  do ele_a = 1, size(enlist_a, 2)
    tet_a%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tet_b%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tets(tet_a, tet_b, tets_c, n_tets_c)
      do ele_c = 1, n_tets_c
        call add(volume_c, tetrahedron_volume(tets_c(ele_c)%v))
      end do
    end do    
  end do
  call report_test("[intersect_tets]", sum(volume_c) .fne. 1.0_real_kind / 3.0_real_kind, .false., "Incorrect intersection volume")

  call initialise(volume_c)
  do ele_a = 1, size(enlist_a, 2)
    tet_a_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_tets(tet_a_real, positions_b(:, enlist_b(:, ele_b)), tets_c_real, n_tets_c)
      do ele_c = 1, n_tets_c
        call add(volume_c, tetrahedron_volume(tets_c_real(:, :, ele_c)))
      end do
    end do    
  end do
  call report_test("[intersect_tets]", sum(volume_c) .fne. 1.0_real_kind / 3.0_real_kind, .false., "Incorrect intersection volume")

  call initialise(volume_c)
  do ele_a = 1, size(enlist_a, 2)
    tet_a%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tet_b%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tets(tet_a, get_planes(tet_b), tets_c, n_tets_c)
      do ele_c = 1, n_tets_c
        call add(volume_c, tetrahedron_volume(tets_c(ele_c)%v))
      end do
    end do    
  end do
  call report_test("[intersect_tets]", sum(volume_c) .fne. 1.0_real_kind / 3.0_real_kind, .false., "Incorrect intersection volume")

  call initialise(volume_c)
  do ele_a = 1, size(enlist_a, 2)
    tet_a%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tet_b%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_polys(tet_a, get_planes(tet_b), tets_c, n_tets_c, work = work)
      do ele_c = 1, n_tets_c
        call add(volume_c, tetrahedron_volume(tets_c(ele_c)%v))
      end do
    end do    
  end do
  call report_test("[intersect_polys]", sum(volume_c) .fne. 1.0_real_kind / 3.0_real_kind, .false., "Incorrect intersection volume")
  
  call initialise(volume_c)
  do ele_a = 1, size(enlist_a, 2)
    tet_a_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_simplices(tet_a_real, positions_b(:, enlist_b(:, ele_b)), tets_c_real, n_tets_c)
      do ele_c = 1, n_tets_c
        call add(volume_c, tetrahedron_volume(tets_c_real(:, :, ele_c)))
      end do
    end do    
  end do
  call report_test("[intersect_simplices]", sum(volume_c) .fne. 1.0_real_kind / 3.0_real_kind, .false., "Incorrect intersection volume")
  
  call initialise(volume_c)
  do ele_a = 1, size(enlist_a, 2)
    tet_a_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(tet_a_real, positions_b(:, enlist_b(:, ele_b)), tets_c_real, n_tets_c)
      do ele_c = 1, n_tets_c
        call add(volume_c, tetrahedron_volume(tets_c_real(:, :, ele_c)))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", sum(volume_c) .fne. 1.0_real_kind / 3.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_tet_intersector
