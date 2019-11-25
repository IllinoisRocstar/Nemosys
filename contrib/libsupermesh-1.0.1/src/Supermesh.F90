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
! femtools/Supermesh.F90 in Fluidity git revision
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

module libsupermesh_supermesh

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_interval_intersection, only : intersect_intervals, &
    & interval_buf_size, interval_size
  use libsupermesh_precision, only : real_kind
  use libsupermesh_tri_intersection, only : tri_type, line_type, max_n_tris_c, &
    & intersect_tris, tri_buf_size, intersect_polys, get_lines, triangle_area
  use libsupermesh_tet_intersection, only : tet_type, plane_type, &
    & max_n_tets_c, intersect_tets, tet_buf_size, intersect_polys, &
    & tetrahedron_volume

  implicit none

  private

  public :: intersect_intervals, interval_buf_size, interval_size
  public :: tri_type, line_type, max_n_tris_c, intersect_tris, tri_buf_size, &
    & intersect_polys, get_lines, triangle_area
  public :: tet_type, plane_type, max_n_tets_c, intersect_tets, tet_buf_size, &
    & tetrahedron_volume
  public :: intersect_tri_quad, intersect_quads, intersect_tet_hex, &
    & intersect_hexes
  public :: max_n_simplices_c, intersect_simplices, intersect_elements, &
    & simplex_volume
  public :: divide_polygon, divide_hex, divide_prism, divide_pyramid
  
  interface max_n_simplices_c
    module procedure max_n_simplices_c_simplices, max_n_simplices_c_elements
  end interface max_n_simplices_c

contains

  pure elemental function max_n_simplices_c_simplices(dim) result(size)
    integer, intent(in) :: dim

    integer :: size

    select case(dim)
      case(1)
        size = interval_buf_size
      case(2)
        size = tri_buf_size
      case(3)
        size = tet_buf_size
      case default
        size = 0
    end select

  end function max_n_simplices_c_simplices

  pure elemental function max_n_simplices_c_elements(dim, loc_a, loc_b) result(size)
    integer, intent(in) :: dim
    integer, intent(in) :: loc_a
    integer, intent(in) :: loc_b

    integer :: size

    if(loc_a == dim + 1 .and. loc_b == dim + 1) then
      size = max_n_simplices_c(dim)
    else if((loc_a == dim + 1 .and. loc_b == 2 ** dim) &
     & .or. (loc_b == dim + 1 .and. loc_a == 2 ** dim)) then
      select case(dim)
        case(2)
          size = max_n_tris_c(n_lines_b = 4)
        case(3)
          size = max_n_tets_c(n_planes_b = 6)
        case default
          size = 0
      end select
    else if(loc_a == 2 ** dim .and. loc_b == 2 ** dim) then
      select case(dim)
        case(2)
          size = max_n_tris_c(n_lines_a = 4, n_lines_b = 4)
        case(3)
          size = 5 * max_n_tets_c(n_planes_b = 6)
        case default
          size = 0
      end select
    else
      size = 0
    end if

  end function max_n_simplices_c_elements

  subroutine intersect_simplices(simplex_a, simplex_b, simplices_c, n_simplices_c)
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: simplex_a
    ! dim x loc_b
    real(kind = real_kind), dimension(:, :), intent(in) :: simplex_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: simplices_c
    integer, intent(out) :: n_simplices_c
    
    select case(size(simplex_a, 1))
      case(1)
        call intersect_intervals(simplex_a(1, :), simplex_b(1, :), simplices_c(1, :, 1), n_simplices_c)
      case(2)
        call intersect_tris(simplex_a, simplex_b, simplices_c, n_simplices_c)
      case(3)
        call intersect_tets(simplex_a, simplex_b, simplices_c, n_simplices_c)
      case default
        libsupermesh_abort("Unsupported element type")
    end select
    
  end subroutine intersect_simplices

  subroutine intersect_elements(element_a, element_b, elements_c, nelements_c)
    ! dim x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: element_a
    ! dim x loc_b
    real(kind = real_kind), dimension(:, :), intent(in) :: element_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: elements_c
    integer, intent(out) :: nelements_c

    integer :: dim, loc_a, loc_b
  
    dim = size(element_a, 1)
    loc_a = size(element_a, 2)
    loc_b = size(element_b, 2)
    
    if(loc_a == dim + 1 .and. loc_b == dim + 1) then
      call intersect_simplices(element_a, element_b, elements_c, nelements_c)
    else if(loc_a == dim + 1 .and. loc_b == 2 ** dim) then
      select case(dim)
        case(2)
          call intersect_tri_quad(element_a, element_b, elements_c, nelements_c)
        case(3)
          call intersect_tet_hex(element_a, element_b, elements_c, nelements_c)
        case default
          libsupermesh_abort("Unsupported element type")
      end select
    else if(loc_a == 2 ** dim .and. loc_b == dim + 1) then
      select case(dim)
        case(2)
          call intersect_tri_quad(element_b, element_a, elements_c, nelements_c)
        case(3)
          call intersect_tet_hex(element_b, element_a, elements_c, nelements_c)
        case default
          libsupermesh_abort("Unsupported element type")
      end select
    else if(loc_a == 2 ** dim .and. loc_b == 2 ** dim) then
      select case(dim)
        case(2)
          call intersect_quads(element_a, element_b, elements_c, nelements_c)
        case(3)
          call intersect_hexes(element_a, element_b, elements_c, nelements_c)
        case default
          libsupermesh_abort("Unsupported element type")
      end select
    else
      libsupermesh_abort("Unsupported element type")
    end if

  end subroutine intersect_elements
  
  subroutine intersect_tri_quad(tri_a, quad_b, tris_c, n_tris_c)
    ! 2 x 3
    real(kind = real_kind), dimension(:, :), intent(in) :: tri_a
    ! 2 x 4
    real(kind = real_kind), dimension(:, :), intent(in) :: quad_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: tris_c
    integer, intent(out) :: n_tris_c
    
    integer :: i
    real(kind = real_kind), dimension(2, 3 * (2 ** 4), 2), save :: work
    type(tri_type) :: ltri_a
    type(tri_type), dimension(3 * (2 ** 4) - 2), save :: ltris_c
    
    ltri_a%v = tri_a
    call intersect_polys(ltri_a, get_lines(quad_b), ltris_c, n_tris_c, work = work)
    do i = 1, n_tris_c
      tris_c(:, :, i) = ltris_c(i)%v
    end do
  
  end subroutine intersect_tri_quad

  subroutine intersect_quads(quad_a, quad_b, tris_c, n_tris_c)
    ! 2 x 4
    real(kind = real_kind), dimension(:, :), intent(in) :: quad_a
    ! 2 x 4
    real(kind = real_kind), dimension(:, :), intent(in) :: quad_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: tris_c
    integer, intent(out) :: n_tris_c
    
    integer :: i
    real(kind = real_kind), dimension(2, 4 * (2 ** 4), 2), save :: work
    type(tri_type), dimension(4 * (2 ** 4) - 2), save :: ltris_c

    call intersect_polys(quad_a, quad_b, ltris_c, n_tris_c, work = work)
    do i = 1, n_tris_c
      tris_c(:, :, i) = ltris_c(i)%v
    end do
  
  end subroutine intersect_quads

  subroutine intersect_tet_hex(tet_a, hex_b, tets_c, n_tets_c)
    ! 3 x 4
    real(kind = real_kind), dimension(:, :), intent(in) :: tet_a
    ! 4 x 8
    real(kind = real_kind), dimension(:, :), intent(in) :: hex_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: tets_c
    integer, intent(out) :: n_tets_c
    
    integer :: i, ln_tets_c
    type(plane_type), dimension(6) :: planes_b
    type(tet_type) :: ltet_a
    type(tet_type), dimension(3 ** 6), save :: ltets_c, work
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    ltet_a%v = tet_a
    planes_b(1)%normal = unit_cross(hex_b(:, 2) - hex_b(:, 1), hex_b(:, 6) - hex_b(:, 1))
    planes_b(1)%c = dot_product(hex_b(:, 1), planes_b(1)%normal)
    planes_b(2)%normal = unit_cross(hex_b(:, 6) - hex_b(:, 5), hex_b(:, 7) - hex_b(:, 5))
    planes_b(2)%c = dot_product(hex_b(:, 5), planes_b(2)%normal)
    planes_b(3)%normal = unit_cross(hex_b(:, 2) - hex_b(:, 6), hex_b(:, 3) - hex_b(:, 6))
    planes_b(3)%c = dot_product(hex_b(:, 6), planes_b(3)%normal)
    planes_b(4)%normal = unit_cross(hex_b(:, 1) - hex_b(:, 2), hex_b(:, 4) - hex_b(:, 2))
    planes_b(4)%c = dot_product(hex_b(:, 2), planes_b(4)%normal)
    planes_b(5)%normal = unit_cross(hex_b(:, 5) - hex_b(:, 1), hex_b(:, 8) - hex_b(:, 1))
    planes_b(5)%c = dot_product(hex_b(:, 1), planes_b(5)%normal)
    planes_b(6)%normal = unit_cross(hex_b(:, 7) - hex_b(:, 8), hex_b(:, 3) - hex_b(:, 8))
    planes_b(6)%c = dot_product(hex_b(:, 8), planes_b(6)%normal)
    
    call intersect_polys(ltet_a, planes_b, ltets_c, ln_tets_c, work = work)
    do i = 1, ln_tets_c
      tets_c(:, :, i) = ltets_c(i)%v
    end do
    n_tets_c = ln_tets_c
   
  contains

    pure function unit_cross(vec_a, vec_b) result(cross)
      real(kind = real_kind), dimension(3), intent(in) :: vec_a, vec_b
      
      real(kind = real_kind), dimension(3) :: cross
      
      cross(1) = vec_a(2) * vec_b(3) - vec_a(3) * vec_b(2)
      cross(2) = vec_a(3) * vec_b(1) - vec_a(1) * vec_b(3)
      cross(3) = vec_a(1) * vec_b(2) - vec_a(2) * vec_b(1)
      cross = cross / norm2(cross)
      
    end function unit_cross
   
  end subroutine intersect_tet_hex

  subroutine intersect_hexes(hex_a, hex_b, tets_c, n_tets_c)
    ! 3 x 8
    real(kind = real_kind), dimension(:, :), intent(in) :: hex_a
    ! 3 x 8
    real(kind = real_kind), dimension(:, :), intent(in) :: hex_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: tets_c
    integer, intent(out) :: n_tets_c
    
    integer :: ln_tets_c
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    ! Cube slicing as e.g. in:
    !   Isosurfaces: Geometry, Topology, and Algorithms, R. Wenger, CRC Press,
    !   2013, figure 2.29
    
    ! Slicing off two opposite corners on the top ...
    call intersect_tet_hex(hex_a(:, (/3, 6, 7, 8/)), hex_b, tets_c, ln_tets_c)
    n_tets_c = ln_tets_c
    call intersect_tet_hex(hex_a(:, (/1, 3, 4, 8/)), hex_b, tets_c(:, :, 1 + n_tets_c:), ln_tets_c)
    n_tets_c = n_tets_c + ln_tets_c
    ! ... and two opposite corners on the bottom ...
    call intersect_tet_hex(hex_a(:, (/1, 5, 6, 8/)), hex_b, tets_c(:, :, 1 + n_tets_c:), ln_tets_c)
    n_tets_c = n_tets_c + ln_tets_c
    call intersect_tet_hex(hex_a(:, (/1, 2, 3, 6/)), hex_b, tets_c(:, :, 1 + n_tets_c:), ln_tets_c)
    n_tets_c = n_tets_c + ln_tets_c
    ! ... to leave a single tetrahedron in the centre
    call intersect_tet_hex(hex_a(:, (/1, 3, 6, 8/)), hex_b, tets_c(:, :, 1 + n_tets_c:), ln_tets_c)
    n_tets_c = n_tets_c + ln_tets_c
  
  end subroutine intersect_hexes

  pure function simplex_volume(simplex) result(volume)
    ! dim x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: simplex

    real(kind = real_kind) :: volume

    select case(size(simplex, 1))
      case(1)
        volume = interval_size(simplex)
      case(2)
        volume = triangle_area(simplex)
      case(3)
        volume = tetrahedron_volume(simplex)
      case default
        volume = -huge(volume)
    end select
  
  end function simplex_volume
  
  pure function divide_polygon(poly) result(tris)
    ! 2 x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: poly
    
    type(tri_type), dimension(size(poly, 2) - 2) :: tris
  
    integer :: i
  
    ! Assumes clockwise or anti-clockwise ordering
    
    forall(i = 1:(size(poly, 2) - 2))
      tris(i)%v(:, 1) = poly(:, 1)
      tris(i)%v(:, 2) = poly(:, i + 1)
      tris(i)%v(:, 3) = poly(:, i + 2)
    end forall
    
  end function divide_polygon
  
  pure function divide_hex(hex) result(tets)
    real(kind = real_kind), dimension(3, 8), intent(in) :: hex
    
    type(tet_type), dimension(5) :: tets
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    ! Cube slicing as e.g. in:
    !   Isosurfaces: Geometry, Topology, and Algorithms, R. Wenger, CRC Press,
    !   2013, figure 2.29
    
    ! Slicing off two opposite corners on the top ...
    tets(1)%v = hex(:, (/3, 6, 7, 8/))
    tets(2)%v = hex(:, (/1, 3, 4, 8/))
    ! ... and two opposite corners on the bottom ...
    tets(3)%v = hex(:, (/1, 5, 6, 8/))
    tets(4)%v = hex(:, (/1, 2, 3, 6/))
    ! ... to leave a single tetrahedron in the centre
    tets(5)%v = hex(:, (/1, 3, 6, 8/))
    
  end function divide_hex
  
  pure function divide_prism(prism) result(tets)
    real(kind = real_kind), dimension(3, 6), intent(in) :: prism
    
    type(tet_type), dimension(3) :: tets
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    tets(1)%v = prism(:, (/1, 2, 3, 5/))
    tets(2)%v = prism(:, (/1, 3, 5, 6/))
    tets(3)%v = prism(:, (/1, 4, 5, 6/))
  
  end function divide_prism
  
  pure function divide_pyramid(pyramid) result(tets)
    real(kind = real_kind), dimension(3, 5), intent(in) :: pyramid
    
    type(tet_type), dimension(2) :: tets
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    tets(1)%v = pyramid(:, (/1, 2, 3, 5/))
    tets(2)%v = pyramid(:, (/1, 3, 4, 5/))
  
  end function divide_pyramid

end module libsupermesh_supermesh
