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
! femtools/tests/test_intersection_finder_2d.F90 in Fluidity git revision
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

subroutine test_intersection_finder_2d() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate
  use libsupermesh_intersection_finder, only : &
    & advancing_front_intersection_finder
  use libsupermesh_precision, only : real_kind
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_unittest_tools, only : report_test

  implicit none

  integer :: i
  logical :: fail
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real(kind = real_kind), dimension(:, :), allocatable :: positions_a, &
    & positions_b
  type(intersections), dimension(1) :: map_ab
  type(intersections), dimension(3) :: bigger_map_ab

  integer, parameter :: dim = 2
  
  call read_node("data/triangle.1.node", dim, positions_a)
  call read_ele("data/triangle.1.ele", dim, enlist_a)
  call read_node("data/triangle.1.node", dim, positions_b)
  call read_ele("data/triangle.1.ele", dim, enlist_b)
  call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  fail = (map_ab(1)%n /= 1)
  call report_test("[intersection finder: length]", fail, .false., "There shall be only one element")

  i = map_ab(1)%v(1)
  fail = (i /= 1)
  call report_test("[intersection finder: correct]", fail, .false., "The answer should be one")

  call deallocate(map_ab)
  deallocate(positions_b, enlist_b)
  call read_node("data/triangle.2.node", dim, positions_b)
  call read_ele("data/triangle.2.ele", dim, enlist_b)
  call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  fail = (map_ab(1)%n /= 3)
  call report_test("[intersection finder: length]", fail, .false., "There shall be three elements")
  do i = 1, 3
    fail = (.not. any(map_ab(1)%v(:map_ab(1)%n) == i))
    call report_test("[intersection finder: correct]", fail, .false., "The answer should be correct")
  end do

  call deallocate(map_ab)
  deallocate(positions_a, enlist_a)
  call read_node("data/triangle.2.node", dim, positions_a)
  call read_ele("data/triangle.2.ele", dim, enlist_a)
  
  call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, bigger_map_ab)
  do i = 1, 3
    fail = (bigger_map_ab(i)%n < 1)
    call report_test("[intersection finder: length]", fail, .false., "There shall be at least one element")

    fail = (.not. any(bigger_map_ab(i)%v(:bigger_map_ab(i)%n) == i))
    call report_test("[intersection finder: correct]", fail, .false., "The answer should be correct")
  end do

  call deallocate(bigger_map_ab)
  deallocate(positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_intersection_finder_2d
