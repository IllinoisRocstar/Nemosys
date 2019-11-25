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

subroutine test_element_intersector() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_precision, only : real_kind
  use libsupermesh_supermesh, only : max_n_simplices_c, intersect_elements, &
    & simplex_volume
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)
  
  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_elements_c
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real(kind = real_kind) :: volume_c
  real(kind = real_kind), dimension(:, :), allocatable :: element_a, &
    & positions_a, positions_b
  real(kind = real_kind), dimension(:, :, :), allocatable :: elements_c
  type(intersections), dimension(:), allocatable :: map_ab
  
  ! Interval-interval
  
  call read_node("data/line.1.node", 1, positions_a)
  call read_ele("data/line.1.ele", 1, enlist_a)
  call read_node("data/line.2.node", 1, positions_b)
  call read_ele("data/line.2.ele", 1, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Tri-tri
  
  call read_node("data/square.1.node", 2, positions_a)
  call read_ele("data/square.1.ele", 2, enlist_a)
  call read_node("data/square.2.node", 2, positions_b)
  call read_ele("data/square.2.ele", 2, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Tri-quad
  
  call read_node("data/square.1.node", 2, positions_a)
  call read_ele("data/square.1.ele", 2, enlist_a)
  call read_node("data/square.4.node", 2, positions_b)
  call read_ele("data/square.4.ele", 2, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Quad-tri
  
  call read_node("data/square.3.node", 2, positions_a)
  call read_ele("data/square.3.ele", 2, enlist_a)
  call read_node("data/square.2.node", 2, positions_b)
  call read_ele("data/square.2.ele", 2, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Quad-quad
  
  call read_node("data/square.3.node", 2, positions_a)
  call read_ele("data/square.3.ele", 2, enlist_a)
  call read_node("data/square.4.node", 2, positions_b)
  call read_ele("data/square.4.ele", 2, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Tet-tet
  
  call read_node("data/cube.1.node", 3, positions_a)
  call read_ele("data/cube.1.ele", 3, enlist_a)
  call read_node("data/cube.2.node", 3, positions_b)
  call read_ele("data/cube.2.ele", 3, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Tet-hex
  
  call read_node("data/cube.1.node", 3, positions_a)
  call read_ele("data/cube.1.ele", 3, enlist_a)
  call read_node("data/cube.3.node", 3, positions_b)
  call read_ele("data/cube.3.ele", 3, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 0.5_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Hex-tet
  
  call read_node("data/cube.3.node", 3, positions_a)
  call read_ele("data/cube.3.ele", 3, enlist_a)
  call read_node("data/cube.2.node", 3, positions_b)
  call read_ele("data/cube.2.ele", 3, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 0.5_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)
  
  ! Hex-hex
  
  call read_node("data/cube.3.node", 3, positions_a)
  call read_ele("data/cube.3.ele", 3, enlist_a)
  call read_node("data/cube.4.node", 3, positions_b)
  call read_ele("data/cube.4.ele", 3, enlist_b)
  allocate(element_a(size(positions_a, 1), size(enlist_a, 1)), &
    & elements_c(size(positions_a, 1), size(positions_a, 1) + 1, &
               & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0_real_kind
  do ele_a = 1, size(enlist_a, 2)
    element_a = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(element_a, positions_b(:, enlist_b(:, ele_b)), elements_c, n_elements_c)
      do ele_c = 1, n_elements_c
        volume_c = volume_c + simplex_volume(elements_c(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 0.5_real_kind, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, element_a, elements_c)

end subroutine test_element_intersector
