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

#define BUF_SIZE 81

module libsupermesh_tet_intersection

  use libsupermesh_precision, only : real_kind

  implicit none

  private

  public :: tet_type, plane_type, max_n_tets_c, intersect_tets, &
    & intersect_polys, get_planes, tetrahedron_volume

  type tet_type
    real(kind = real_kind), dimension(3, 4) :: v ! vertices of the tet
    integer, dimension(4) :: colours = -1 ! surface colours
  end type tet_type

  type plane_type
    real(kind = real_kind), dimension(3) :: normal
    real(kind = real_kind) :: c
  end type plane_type

  interface max_n_tets_c
    module procedure max_n_tets_c_tet
  end interface max_n_tets_c

  interface intersect_tets
    module procedure intersect_tets_real, intersect_tets_tet, &
      & intersect_tets_planes
  end interface intersect_tets

  interface intersect_polys
    module procedure intersect_tet_planes
  end interface intersect_polys

  interface get_planes
    module procedure get_planes_tet
  end interface get_planes
  
  interface tetrahedron_volume
    module procedure tetrahedron_volume_real, tetrahedron_volume_tet
  end interface tetrahedron_volume

  integer, parameter, public :: tet_buf_size = BUF_SIZE

  type(tet_type), dimension(BUF_SIZE), save :: tets_tmp_buf
  integer, save :: n_tets_tmp = 0

contains

  subroutine intersect_tets_real(tet_a, tet_b, tets_c, n_tets_c)
    real(kind = real_kind), dimension(3, 4), intent(in) :: tet_a
    real(kind = real_kind), dimension(3, 4), intent(in) :: tet_b
    real(kind = real_kind), dimension(:, :, :), intent(inout) :: tets_c
    integer, intent(out) :: n_tets_c

    integer :: i
    type(tet_type) :: tet_a_t, tet_b_t
    type(tet_type), dimension(BUF_SIZE), save :: tets_c_t

    tet_a_t%v = tet_a
    tet_b_t%v = tet_b
    call intersect_tets(tet_a_t, tet_b_t, tets_c_t, n_tets_c)

    do i = 1, n_tets_c
      tets_c(:, :, i) = tets_c_t(i)%v
    end do

  end subroutine intersect_tets_real

  subroutine intersect_tets_tet(tet_a, tet_b, tets_c, n_tets_c)
    type(tet_type), intent(in) :: tet_a
    type(tet_type), intent(in) :: tet_b
    type(tet_type), dimension(BUF_SIZE), intent(inout) :: tets_c
    integer, intent(out) :: n_tets_c

    call intersect_tets_planes(tet_a, get_planes(tet_b), tets_c, n_tets_c, vol_b = tetrahedron_volume(tet_b))

  end subroutine intersect_tets_tet

  subroutine intersect_tets_planes(tet_a, planes_b, tets_c, n_tets_c, vol_b)
    type(tet_type), intent(in) :: tet_a
    type(plane_type), dimension(4), intent(in)  :: planes_b
    type(tet_type), dimension(BUF_SIZE), intent(inout) :: tets_c
    integer, intent(out) :: n_tets_c
    real(kind = real_kind), optional, intent(in) :: vol_b

    integer :: i, j
    real(kind = real_kind) :: tol, vol

    n_tets_c = 1
    tets_c(1) = tet_a

    if(present(vol_b)) then
      tol = 10.0_real_kind * min(spacing(tetrahedron_volume(tet_a)), spacing(vol_b))
    else
      tol = 10.0_real_kind * spacing(tetrahedron_volume(tet_a))
    end if
    do i = 1, size(planes_b)
      ! Clip the tet_array against the i'th plane
      n_tets_tmp = 0

      do j = 1, n_tets_c
        call clip_buf(planes_b(i), tets_c(j))
      end do

      if(i /= size(planes_b)) then
        n_tets_c = n_tets_tmp
        tets_c(:n_tets_c) = tets_tmp_buf(:n_tets_c)
      else
        ! Copy the result if the volume is >= tol
        n_tets_c = 0
        do j = 1, n_tets_tmp
          vol = tetrahedron_volume(tets_tmp_buf(j))
          if(vol >= tol) then
            n_tets_c = n_tets_c + 1
            tets_c(n_tets_c) = tets_tmp_buf(j)
          end if
        end do
      end if
    end do

  end subroutine intersect_tets_planes
  
  pure elemental function max_n_tets_c_tet(n_planes_b) result(max_n_tets_c)
    integer, intent(in) :: n_planes_b
    
    integer :: max_n_tets_c
    
    max_n_tets_c = 3 ** n_planes_b
    
  end function max_n_tets_c_tet
  
  pure subroutine intersect_tet_planes(tet_a, planes_b, tets_c, n_tets_c, vol_b, work)
    type(tet_type), intent(in) :: tet_a
    type(plane_type), dimension(:), intent(in)  :: planes_b
    type(tet_type), dimension(:), intent(inout) :: tets_c
    integer, intent(out) :: n_tets_c
    real(kind = real_kind), optional, intent(in) :: vol_b
    type(tet_type), dimension(:), target, optional, intent(inout) :: work

    integer :: i, j, ntets_new
    real(kind = real_kind) :: tol, vol
    type(tet_type), dimension(:), pointer :: tets_new
    
    if(present(work)) then
      tets_new => work
    else
      allocate(tets_new(max_n_tets_c(size(planes_b))))
    end if

    n_tets_c = 1
    tets_c(1) = tet_a

    if(present(vol_b)) then
      tol = 10.0_real_kind * min(spacing(tetrahedron_volume(tet_a)), spacing(vol_b))
    else
      tol = 10.0_real_kind * spacing(tetrahedron_volume(tet_a))
    end if
    do i = 1, size(planes_b)
      ! Clip the tet_array against the i'th plane
      ntets_new = 0

      do j = 1, n_tets_c
        call clip(planes_b(i), tets_c(j), tets_new, ntets_new)
      end do

      if(i /= size(planes_b)) then
        n_tets_c = ntets_new
        tets_c(:n_tets_c) = tets_new(:n_tets_c)
      else
        ! Copy the result if the volume is >= tol
        n_tets_c = 0
        do j = 1, ntets_new
          vol = tetrahedron_volume(tets_new(j))
          if(vol >= tol) then
            n_tets_c = n_tets_c + 1
            tets_c(n_tets_c) = tets_new(j)
          end if
        end do
      end if
    end do
    
    if(.not. present(work)) deallocate(tets_new)
  
  end subroutine intersect_tet_planes

  subroutine clip_buf(plane, tet)
    ! Clip tet against the plane and append any output to tets_tmp_buf.
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet

    real(kind = real_kind), dimension(4) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(4) :: neg_idx, pos_idx, zer_idx
    integer :: i

    real(kind = real_kind) :: invdiff, w0, w1
    type(tet_type) :: tet_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    zer_cnt = 0

    dists = distances_to_plane(plane, tet)
    do i=1,4
      if (abs(dists(i)) < epsilon(dists(i))) then
        zer_cnt = zer_cnt + 1
        zer_idx(zer_cnt) = i
      else if (dists(i) < 0.0_real_kind) then
        neg_cnt = neg_cnt + 1
        neg_idx(neg_cnt) = i
      else if (dists(i) > 0.0_real_kind) then
        pos_cnt = pos_cnt + 1
        pos_idx(pos_cnt) = i
      end if
    end do

    if (neg_cnt == 0) then
      ! tet is completely on positive side of plane, full clip
      return
    end if

    if (pos_cnt == 0) then
      ! tet is completely on negative side of plane, no clip
      n_tets_tmp = n_tets_tmp + 1
      tets_tmp_buf(n_tets_tmp) = tet
      return
    end if

    ! The tet is split by the plane, so we have more work to do.

    select case(pos_cnt)
    case(3)
      ! +++-
      n_tets_tmp = n_tets_tmp + 1
      tets_tmp_buf(n_tets_tmp) = tet
      do i=1,pos_cnt
        invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(i)) * invdiff
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) = &
           w0 * tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) + &
           w1 * tets_tmp_buf(n_tets_tmp)%v(:, neg_idx(1))
      end do
      ! The colours will have been inherited already; we just need to zero
      ! the one corresponding to the plane cut
      tets_tmp_buf(n_tets_tmp)%colours(face_no(pos_idx(1), pos_idx(2), pos_idx(3))) = 0
    case(2)
      select case(neg_cnt)
      case(2)
        ! ++--
        do i=1,pos_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(1))
        end do
        do i=1,neg_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(2)) )
          w0 = -dists(neg_idx(2)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i+2) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(2))
        end do

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(2)) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(2)) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet%v(:, neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet_tmp%v(:, 4)
        tets_tmp_buf(n_tets_tmp)%colours(2) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%colours(3) = tet%colours(pos_idx(1))
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(1))

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet%v(:, neg_idx(1))
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(2) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(3) = tet%colours(pos_idx(2))
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(2))
      case(1)
        ! ++-0
        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        do i=1,pos_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) = &
             w0 * tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) + &
             w1 * tets_tmp_buf(n_tets_tmp)%v(:, neg_idx(1))
        end do
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0
      end select
    case(1)
      select case(neg_cnt)
      case(3)
        ! +---
        do i=1,neg_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(1) = tet%colours(neg_idx(1))
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet%v(:, neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%colours(2) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet%v(:, neg_idx(3))
        tets_tmp_buf(n_tets_tmp)%colours(3) = tet%colours(neg_idx(3))
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(4) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet%v(:, neg_idx(3))
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(2) = tet%colours(neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%colours(3) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(1))
      case(2)
        ! +--0
        do i=1,neg_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet%v(:, zer_idx(1))
        tets_tmp_buf(n_tets_tmp)%colours(2) = tet%colours(zer_idx(1))
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet%v(:, neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%colours(3) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(1))
      case(1)
        ! +-00
        invdiff = 1.0_real_kind / ( dists(pos_idx(1)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(1)) * invdiff

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(1))
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0
      end select
    end select

  end subroutine clip_buf

  pure subroutine clip(plane, tet, tets_new, ntets_new)
    ! Clip tet against the plane.
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet
    type(tet_type), dimension(:), intent(inout) :: tets_new
    integer, intent(inout) :: ntets_new

    real(kind = real_kind), dimension(4) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(4) :: neg_idx, pos_idx, zer_idx
    integer :: i

    real(kind = real_kind) :: invdiff, w0, w1
    type(tet_type) :: tet_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    zer_cnt = 0

    dists = distances_to_plane(plane, tet)
    do i=1,4
      if (abs(dists(i)) < epsilon(dists(i))) then
        zer_cnt = zer_cnt + 1
        zer_idx(zer_cnt) = i
      else if (dists(i) < 0.0_real_kind) then
        neg_cnt = neg_cnt + 1
        neg_idx(neg_cnt) = i
      else if (dists(i) > 0.0_real_kind) then
        pos_cnt = pos_cnt + 1
        pos_idx(pos_cnt) = i
      end if
    end do

    if (neg_cnt == 0) then
      ! tet is completely on positive side of plane, full clip
      return
    end if

    if (pos_cnt == 0) then
      ! tet is completely on negative side of plane, no clip
      ntets_new = ntets_new + 1
      tets_new(ntets_new) = tet
      return
    end if

    ! The tet is split by the plane, so we have more work to do.

    select case(pos_cnt)
    case(3)
      ! +++-
      ntets_new = ntets_new + 1
      tets_new(ntets_new) = tet
      do i=1,pos_cnt
        invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(i)) * invdiff
        tets_new(ntets_new)%v(:, pos_idx(i)) = &
           w0 * tets_new(ntets_new)%v(:, pos_idx(i)) + &
           w1 * tets_new(ntets_new)%v(:, neg_idx(1))
      end do
      ! The colours will have been inherited already; we just need to zero
      ! the one corresponding to the plane cut
      tets_new(ntets_new)%colours(face_no(pos_idx(1), pos_idx(2), pos_idx(3))) = 0
    case(2)
      select case(neg_cnt)
      case(2)
        ! ++--
        do i=1,pos_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(1))
        end do
        do i=1,neg_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(2)) )
          w0 = -dists(neg_idx(2)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i+2) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(2))
        end do

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%v(:, pos_idx(2)) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(neg_idx(1)) = 0
        tets_new(ntets_new)%colours(neg_idx(2)) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet%v(:, neg_idx(2))
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet_tmp%v(:, 4)
        tets_new(ntets_new)%colours(2) = 0
        tets_new(ntets_new)%v(:, 3) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%colours(3) = tet%colours(pos_idx(1))
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(1))

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet%v(:, neg_idx(1))
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(2) = 0
        tets_new(ntets_new)%v(:, 3) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(3) = tet%colours(pos_idx(2))
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(2))
      case(1)
        ! ++-0
        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        do i=1,pos_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tets_new(ntets_new)%v(:, pos_idx(i)) = &
             w0 * tets_new(ntets_new)%v(:, pos_idx(i)) + &
             w1 * tets_new(ntets_new)%v(:, neg_idx(1))
        end do
        tets_new(ntets_new)%colours(neg_idx(1)) = 0
      end select
    case(1)
      select case(neg_cnt)
      case(3)
        ! +---
        do i=1,neg_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(neg_idx(1)) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(1) = tet%colours(neg_idx(1))
        tets_new(ntets_new)%v(:, 2) = tet%v(:, neg_idx(2))
        tets_new(ntets_new)%colours(2) = 0
        tets_new(ntets_new)%v(:, 3) = tet%v(:, neg_idx(3))
        tets_new(ntets_new)%colours(3) = tet%colours(neg_idx(3))
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(4) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet%v(:, neg_idx(3))
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(2) = tet%colours(neg_idx(2))
        tets_new(ntets_new)%v(:, 3) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%colours(3) = 0
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(1))
      case(2)
        ! +--0
        do i=1,neg_cnt
          invdiff = 1.0_real_kind / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(neg_idx(1)) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet%v(:, zer_idx(1))
        tets_new(ntets_new)%colours(2) = tet%colours(zer_idx(1))
        tets_new(ntets_new)%v(:, 3) = tet%v(:, neg_idx(2))
        tets_new(ntets_new)%colours(3) = 0
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(1))
      case(1)
        ! +-00
        invdiff = 1.0_real_kind / ( dists(pos_idx(1)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(1)) * invdiff

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(1))
        tets_new(ntets_new)%colours(neg_idx(1)) = 0
      end select
    end select

  end subroutine clip
  
  pure function get_planes_tet(tet) result(planes)
    type(tet_type), intent(in) :: tet
    type(plane_type), dimension(4) :: planes

    real(kind = real_kind), dimension(3) :: edge10, edge20, edge30, edge21, &
      & edge31
    real(kind = real_kind) :: det
    integer :: i

    edge10 = tet%v(:, 2) - tet%v(:, 1);
    edge20 = tet%v(:, 3) - tet%v(:, 1);
    edge30 = tet%v(:, 4) - tet%v(:, 1);
    edge21 = tet%v(:, 3) - tet%v(:, 2);
    edge31 = tet%v(:, 4) - tet%v(:, 2);

    planes(1)%normal = unit_cross(edge20, edge10)
    planes(2)%normal = unit_cross(edge10, edge30)
    planes(3)%normal = unit_cross(edge30, edge20)
    planes(4)%normal = unit_cross(edge21, edge31)

    det = dot_product(edge10, planes(4)%normal)
    if (det < 0) then
      do i=1,4
        planes(i)%normal = -planes(i)%normal
      end do
    end if

    ! And calibrate what is the zero of this plane by dotting with
    ! a point we know to be on it
    do i=1,4
      planes(i)%c = dot_product(tet%v(:, i), planes(i)%normal)
    end do

  end function get_planes_tet

  pure function unit_cross(vec_a, vec_b) result(cross)
    real(kind = real_kind), dimension(3), intent(in) :: vec_a, vec_b
    real(kind = real_kind), dimension(3) :: cross
    cross(1) = vec_a(2) * vec_b(3) - vec_a(3) * vec_b(2)
    cross(2) = vec_a(3) * vec_b(1) - vec_a(1) * vec_b(3)
    cross(3) = vec_a(1) * vec_b(2) - vec_a(2) * vec_b(1)

    cross = cross / norm2(cross)
  end function unit_cross

  pure function distances_to_plane(plane, tet) result(dists)
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet
    real(kind = real_kind), dimension(4) :: dists
    integer :: i

    forall(i=1:4)
      dists(i) = dot_product(plane%normal, tet%v(:, i)) - plane%c
    end forall
  end function distances_to_plane

  pure function tetrahedron_volume_real(tet) result(volume)
    real(kind = real_kind), dimension(3, 4), intent(in) :: tet

    real(kind = real_kind) :: volume

    real(kind = real_kind), dimension(3) :: e1, e2, e3

    e1 = tet(:, 2) - tet(:, 1)
    e2 = tet(:, 3) - tet(:, 1)
    e3 = tet(:, 4) - tet(:, 1)

    volume = (1.0_real_kind / 6.0_real_kind) * &
      & abs(e1(1) * (e2(2) * e3(3) - e2(3) * e3(2)) &
        & + e1(2) * (e2(3) * e3(1) - e2(1) * e3(3)) &
        & + e1(3) * (e2(1) * e3(2) - e2(2) * e3(1)))

  end function tetrahedron_volume_real

  pure elemental function tetrahedron_volume_tet(tet) result(volume)
    type(tet_type), intent(in) :: tet

    real(kind = real_kind) :: volume

    real(kind = real_kind), dimension(3) :: cross, vec_a, vec_b, vec_c

    vec_a = tet%v(:, 1) - tet%v(:, 4)
    vec_b = tet%v(:, 2) - tet%v(:, 4)
    vec_c = tet%v(:, 3) - tet%v(:, 4)

    cross(1) = vec_b(2) * vec_c(3) - vec_b(3) * vec_c(2)
    cross(2) = vec_b(3) * vec_c(1) - vec_b(1) * vec_c(3)
    cross(3) = vec_b(1) * vec_c(2) - vec_b(2) * vec_c(1)

    volume = abs(dot_product(vec_a, cross)) / 6.0_real_kind
    
  end function tetrahedron_volume_tet

  pure elemental function face_no(i, j, k) result(face)
    ! Given three local node numbers, what is the face that they share?
    integer, intent(in) :: i, j, k
    integer :: face

    do face=1,4
      if (face /= i .and. face /= j .and. face /= k) return
    end do

  end function face_no

end module libsupermesh_tet_intersection
