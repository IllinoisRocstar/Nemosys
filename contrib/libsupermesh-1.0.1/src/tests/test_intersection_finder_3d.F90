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
! femtools/tests/test_intersection_finder_3d.F90 in Fluidity git revision
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

subroutine test_intersection_finder_3d() bind(c)

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

  integer, parameter :: dim = 3

  call read_node("data/tet.node", dim, positions_a)
  call read_ele("data/tet.ele", dim, enlist_a)
  call read_node("data/tet.node", dim, positions_b)
  call read_ele("data/tet.ele", dim, enlist_b)
  call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  fail = (map_ab(1)%n /= 1)
  call report_test("[intersection finder: length]", fail, .false., "There shall be only one element")

  i = map_ab(1)%v(1)
  fail = (i /= 1)
  call report_test("[intersection finder: correct]", fail, .false., "The answer should be one")
    
  call deallocate(map_ab)
  deallocate(positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_intersection_finder_3d
