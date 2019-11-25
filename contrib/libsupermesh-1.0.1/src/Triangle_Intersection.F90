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
! femtools/Tetrahedron_intersection.F90 in Fluidity git revision
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

#define BUF_SIZE 22

module libsupermesh_tri_intersection

  use libsupermesh_precision, only : real_kind

  implicit none

  private

  public :: tri_type, line_type, max_n_tris_c, intersect_tris, &
    & intersect_polys, get_lines, triangle_area

  type tri_type
    real(kind = real_kind), dimension(2, 3) :: v
  end type tri_type

  type line_type
    real(kind = real_kind), dimension(2) :: normal
    real(kind = real_kind), dimension(2) :: point
  end type line_type
  
  interface max_n_tris_c
    module procedure max_n_tris_c_tri, max_n_tris_c_poly
  end interface max_n_tris_c

  interface intersect_tris
    module procedure intersect_tris_real, intersect_tris_tri
  end interface intersect_tris
  
  interface intersect_polys
    module procedure intersect_tri_lines, intersect_polys_real, &
      & intersect_polys_lines
  end interface intersect_polys

  interface get_lines
    module procedure get_lines_tri, get_lines_poly
  end interface get_lines
  
  interface triangle_area
    module procedure triangle_area_real, triangle_area_tri
  end interface triangle_area

  integer, parameter, public :: tri_buf_size = BUF_SIZE

  real(kind = real_kind), dimension(2, BUF_SIZE + 2), save :: points_tmp
  integer, save :: n_points_tmp

contains

  pure elemental function max_n_tris_c_tri(n_lines_b) result(max_n_tris_c)
    integer, intent(in) :: n_lines_b
    
    integer :: max_n_tris_c
    
    max_n_tris_c = 3 * (2 ** n_lines_b) - 2
    
  end function max_n_tris_c_tri

  pure elemental function max_n_tris_c_poly(n_lines_a, n_lines_b) result(max_n_tris_c)
    integer, intent(in) :: n_lines_a
    integer, intent(in) :: n_lines_b
    
    integer :: max_n_tris_c
    
    max_n_tris_c = n_lines_a * (2 ** n_lines_b) - 2
    
  end function max_n_tris_c_poly

  subroutine intersect_tris_real(tri_a, tri_b, tris_c, n_tris_c)
    real(kind = real_kind), dimension(2, 3), intent(in) :: tri_a
    real(kind = real_kind), dimension(2, 3), intent(in) :: tri_b
    real(kind = real_kind), dimension(2, 3, BUF_SIZE), intent(out) :: tris_c
    integer, intent(out) :: n_tris_c

    integer :: i
    type(tri_type) :: tri_a_t, tri_b_t
    type(tri_type), dimension(BUF_SIZE), save :: tris_c_t

    tri_a_t%v = tri_a
    tri_b_t%v = tri_b
    call intersect_tris(tri_a_t, tri_b_t, tris_c_t, n_tris_c)
    do i = 1, n_tris_c
      tris_c(:, :, i) = tris_c_t(i)%v
    end do

  end subroutine intersect_tris_real

  subroutine intersect_tris_tri(tri_a, tri_b, tris_c, n_tris_c)
    type(tri_type), intent(in) :: tri_a
    type(tri_type), intent(in) :: tri_b
    type(tri_type), dimension(BUF_SIZE), intent(out) :: tris_c
    integer, intent(out) :: n_tris_c

    integer :: i
    real(kind = real_kind) :: tol
    type(line_type), dimension(3) :: lines_b

    real(kind = real_kind), dimension(2, BUF_SIZE + 2), save :: points
    integer :: n_points

    lines_b = get_lines(tri_b)

    points(:, :3) = tri_a%v
    n_points = 3

    n_tris_c = 0
    call clip_buf(lines_b(1), points, n_points)
    if(n_points_tmp < 3) return
    points(:, :n_points_tmp) = points_tmp(:, :n_points_tmp)
    n_points = n_points_tmp
    call clip_buf(lines_b(2), points, n_points)
    if(n_points_tmp < 3) return
    points(:, :n_points_tmp) = points_tmp(:, :n_points_tmp)
    n_points = n_points_tmp
    call clip_buf(lines_b(3), points, n_points)
    if(n_points_tmp < 3) return

    tol = 10.0_real_kind * min(spacing(triangle_area(tri_a)), spacing(triangle_area(tri_b)))
    do i = 1, n_points_tmp - 2
      n_tris_c = n_tris_c + 1
      tris_c(n_tris_c)%v(:, 1) = points_tmp(:, 1)
      tris_c(n_tris_c)%v(:, 2) = points_tmp(:, i + 1)
      tris_c(n_tris_c)%v(:, 3) = points_tmp(:, i + 2)
      if(triangle_area(tris_c(n_tris_c)) < tol) then
        n_tris_c = n_tris_c - 1
      end if
    end do

  end subroutine intersect_tris_tri

  pure subroutine intersect_tri_lines(tri_a, lines_b, tris_c, n_tris_c, area_b, work)
    type(tri_type), intent(in) :: tri_a
    type(line_type), dimension(:), intent(in) :: lines_b
    type(tri_type), dimension(:), intent(out) :: tris_c
    integer, intent(out) :: n_tris_c
    real(kind = real_kind), optional, intent(in) :: area_b
    real(kind = real_kind), dimension(:, :, :), target, optional, intent(out) :: &
      & work

    call intersect_polys(tri_a%v, lines_b, tris_c, n_tris_c, area_a = triangle_area(tri_a), area_b = area_b, work = work)

  end subroutine intersect_tri_lines
  
  pure subroutine intersect_polys_real(poly_a, poly_b, tris_c, n_tris_c, area_a, area_b, work)
    ! 2 x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: poly_a
    ! 2 x loc_b
    real(kind = real_kind), dimension(:, :), intent(in) :: poly_b
    type(tri_type), dimension(:), intent(out) :: tris_c
    integer, intent(out) :: n_tris_c
    real(kind = real_kind), optional, intent(in) :: area_a
    real(kind = real_kind), optional, intent(in) :: area_b
    real(kind = real_kind), dimension(:, :, :), target, optional, intent(out) :: &
      & work
    
    call intersect_polys(poly_a, get_lines(poly_b), tris_c, n_tris_c, area_a = area_a, area_b = area_b, work = work)
  
  end subroutine intersect_polys_real

  pure subroutine intersect_polys_lines(poly_a, lines_b, tris_c, n_tris_c, area_a, area_b, work)
    ! 2 x loc_a
    real(kind = real_kind), dimension(:, :), intent(in) :: poly_a
    ! loc_b
    type(line_type), dimension(:), intent(in) :: lines_b
    type(tri_type), dimension(:), intent(out) :: tris_c
    integer, intent(out) :: n_tris_c
    real(kind = real_kind), optional, intent(in) :: area_a
    real(kind = real_kind), optional, intent(in) :: area_b
    real(kind = real_kind), dimension(:, :, :), target, optional, intent(out) :: &
      & work

    integer :: i
    real(kind = real_kind) :: tol

    real(kind = real_kind), dimension(:, :), pointer :: points, points_new
    integer :: n_points, n_points_new
    
    if(present(work)) then
      points => work(:, :, 1)
      points_new => work(:, :, 2)
    else
      allocate(points(2, max_n_tris_c(size(poly_a, 2), size(lines_b)) + 2))
      allocate(points_new(2, size(points)))
    end if

    points(:, :size(poly_a, 2)) = poly_a
    n_points = size(poly_a, 2)

    n_tris_c = 0
    do i = 1, size(lines_b)
      call clip(lines_b(i), points, n_points, points_new, n_points_new)
      if(n_points_new < 3) goto 42
      points(:, :n_points_new) = points_new(:, :n_points_new)
      n_points = n_points_new
    end do

    if(present(area_b)) then
      if(present(area_a)) then
        tol = 10.0_real_kind * min(spacing(area_a), spacing(area_b))
      else
        tol = 10.0_real_kind * spacing(area_b)
      end if
    else if(present(area_a)) then
      tol = 10.0_real_kind * spacing(area_a)
    else
      tol = 10.0_real_kind * tiny(tol)
    end if
    do i = 1, n_points_new - 2
      n_tris_c = n_tris_c + 1
      tris_c(n_tris_c)%v(:, 1) = points_new(:, 1)
      tris_c(n_tris_c)%v(:, 2) = points_new(:, i + 1)
      tris_c(n_tris_c)%v(:, 3) = points_new(:, i + 2)
      if(triangle_area(tris_c(n_tris_c)) < tol) then
        n_tris_c = n_tris_c - 1
      end if
    end do

42  if(.not. present(work)) deallocate(points, points_new)

  end subroutine intersect_polys_lines

  ! Sutherland-Hodgman clipping algorithm. See:
  !   Reentrant polygon clipping, I. E. Sutherland and G. W. Hodgman,
  !   Communications of the ACM vol. 17, 1974, pp. 32--42.
  ! Note that the role of the "initial" and "terminal" is swapped here
  subroutine clip_buf(line, points, n_points)
    type(line_type), intent(in) :: line
    real(kind = real_kind), dimension(2, BUF_SIZE + 2), intent(in) :: points
    integer, intent(in) :: n_points

    integer :: i
    real(kind = real_kind) :: d1, d2, f
    real(kind = real_kind), dimension(2) :: p1, p2
    real(kind = real_kind), dimension(BUF_SIZE + 2), save :: d

    do i = 1, n_points
      d(i) = dot_product(line%normal, points(:, i) - line%point)
    end do

    n_points_tmp = 0
    do i = 1, n_points
      p1 = points(:, i)
      d1 = d(i)
      if(i == n_points) then
        p2 = points(:, 1)
        d2 = d(1)
      else
        p2 = points(:, i + 1)
        d2 = d(i + 1)
      end if

      if(d1 < 0.0_real_kind) then
        if(d2 <= 0.0_real_kind) then
          ! No clip
          n_points_tmp = n_points_tmp + 1
          points_tmp(:, n_points_tmp) = p1
        else
          ! New point
          n_points_tmp = n_points_tmp + 1
          points_tmp(:, n_points_tmp) = p1
          n_points_tmp = n_points_tmp + 1
          f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0_real_kind), 0.0_real_kind)
          points_tmp(:, n_points_tmp) = p1 + f * (p2 - p1)
        end if
      else if(d1 == 0.0_real_kind) then
        ! No clip
        n_points_tmp = n_points_tmp + 1
        points_tmp(:, n_points_tmp) = p1
      else if(d2 < 0.0_real_kind) then
        ! Move point
        n_points_tmp = n_points_tmp + 1
        f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0_real_kind), 0.0_real_kind)
        points_tmp(:, n_points_tmp) = p1 + f * (p2 - p1)
      !else
        ! Full clip
      end if
    end do

  end subroutine clip_buf

  ! Sutherland-Hodgman clipping algorithm. See:
  !   Reentrant polygon clipping, I. E. Sutherland and G. W. Hodgman,
  !   Communications of the ACM vol. 17, 1974, pp. 32--42.
  ! Note that the role of the "initial" and "terminal" is swapped here
  pure subroutine clip(line, points, n_points, points_new, n_points_new)
    type(line_type), intent(in) :: line
    real(kind = real_kind), dimension(:, :), intent(in) :: points
    integer, intent(in) :: n_points
    real(kind = real_kind), dimension(:, :), intent(out) :: points_new
    integer, intent(out) :: n_points_new

    integer :: i
    real(kind = real_kind) :: d0, d1, d2, f
    real(kind = real_kind), dimension(2) :: p1, p2

    d0 = dot_product(line%normal, points(:, 1) - line%point)
    d2 = d0
    n_points_new = 0
    do i = 1, n_points
      p1 = points(:, i)
      d1 = d2
      if(i == n_points) then
        p2 = points(:, 1)
        d2 = d0
      else
        p2 = points(:, i + 1)
        d2 = dot_product(line%normal, points(:, i + 1) - line%point)
      end if

      if(d1 < 0.0_real_kind) then
        if(d2 <= 0.0_real_kind) then
          ! No clip
          n_points_new = n_points_new + 1
          points_new(:, n_points_new) = p1
        else
          ! New point
          n_points_new = n_points_new + 1
          points_new(:, n_points_new) = p1
          n_points_new = n_points_new + 1
          f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0_real_kind), 0.0_real_kind)
          points_new(:, n_points_new) = p1 + f * (p2 - p1)
        end if
      else if(d1 == 0.0_real_kind) then
        ! No clip
        n_points_new = n_points_new + 1
        points_new(:, n_points_new) = p1
      else if(d2 < 0.0_real_kind) then
        ! Move point
        n_points_new = n_points_new + 1
        f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0_real_kind), 0.0_real_kind)
        points_new(:, n_points_new) = p1 + f * (p2 - p1)
      !else
        ! Full clip
      end if
    end do

  end subroutine clip

  pure function get_lines_tri(tri) result(lines)
    type(tri_type), intent(in) :: tri
    
    type(line_type), dimension(3) :: lines

    ! Note that the normals are not normalised

    lines(1)%normal(1) = -(tri%v(2, 2) - tri%v(2, 1))
    lines(1)%normal(2) =  (tri%v(1, 2) - tri%v(1, 1))
    lines(1)%point = tri%v(:, 1)

    lines(2)%normal(1) = -(tri%v(2, 3) - tri%v(2, 2))
    lines(2)%normal(2) =  (tri%v(1, 3) - tri%v(1, 2))
    lines(2)%point = tri%v(:, 2)

    lines(3)%normal(1) = -(tri%v(2, 1) - tri%v(2, 3))
    lines(3)%normal(2) =  (tri%v(1, 1) - tri%v(1, 3))
    lines(3)%point = tri%v(:, 3)

    if(dot_product(tri%v(:, 2) - tri%v(:, 1), lines(3)%normal) > 0.0_real_kind) then
      lines(1)%normal = -lines(1)%normal
      lines(2)%normal = -lines(2)%normal
      lines(3)%normal = -lines(3)%normal
    end if

  end function get_lines_tri

  pure function get_lines_poly(poly) result(lines)
    ! 2 x loc
    real(kind = real_kind), dimension(:, :), intent(in) :: poly
    
    type(line_type), dimension(size(poly, 2)) :: lines

    integer :: i, loc
    
    loc = size(poly, 2)
  
    ! Assumes clockwise or anti-clockwise ordering, and that the first three
    ! nodes are not co-linear
    ! Note that the normals are not normalised

    do i = 1, loc - 1
      lines(i)%normal(1) =  (poly(2, i + 1) - poly(2, i))
      lines(i)%normal(2) = -(poly(1, i + 1) - poly(1, i))
      lines(i)%point = poly(:, i)
    end do
    lines(loc)%normal(1) =  (poly(2, 1) - poly(2, loc))
    lines(loc)%normal(2) = -(poly(1, 1) - poly(1, loc))
    lines(loc)%point = poly(:, loc)
    
    if(dot_product(poly(:, 1) - lines(2)%point, lines(2)%normal) > 0.0_real_kind) then
      do i = 1, loc
        lines(i)%normal = -lines(i)%normal
      end do
    end if
 
  end function get_lines_poly

  pure function triangle_area_real(tri) result(area)
    real(kind = real_kind), dimension(2, 3), intent(in) :: tri

    real(kind = real_kind) :: area
    real(kind = real_kind), dimension(2) :: u, v

    u = tri(:, 3) - tri(:, 1)
    v = tri(:, 2) - tri(:, 1)

    area = 0.5_real_kind * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area_real

  pure elemental function triangle_area_tri(tri) result(area)
    type(tri_type), intent(in) :: tri

    real(kind = real_kind) :: area
    real(kind = real_kind), dimension(2) :: u, v

    u = tri%v(:, 3) - tri%v(:, 1)
    v = tri%v(:, 2) - tri%v(:, 1)

    area = 0.5_real_kind * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area_tri
  
end module libsupermesh_tri_intersection
