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
! femtools/Integer_hash_table.F90 in Fluidity git revision
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

module libsupermesh_integer_map

  use iso_c_binding, only : c_int, c_ptr
  
  implicit none

  private
  
  type integer_map
    type(c_ptr), pointer :: ptr
  end type integer_map

  interface
    subroutine cinteger_map_new(i) bind(c, name = "libsupermesh_integer_map_new")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine cinteger_map_new

    subroutine cinteger_map_delete(i) bind(c, name = "libsupermesh_integer_map_delete")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine cinteger_map_delete

    subroutine cinteger_map_insert(i, key, value) bind(c, name = "libsupermesh_integer_map_insert")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
      integer(kind = c_int), value :: value
    end subroutine cinteger_map_insert

    subroutine cinteger_map_size(i, size) bind(c, name = "libsupermesh_integer_map_size")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: size
    end subroutine cinteger_map_size

    subroutine cinteger_map_fetch(i, key, value) bind(c, name = "libsupermesh_integer_map_fetch")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
      integer(kind = c_int) :: value
    end subroutine cinteger_map_fetch

    subroutine cinteger_map_remove(i, key) bind(c, name = "libsupermesh_integer_map_remove")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
    end subroutine cinteger_map_remove

    subroutine cinteger_map_has_key(i, key, present) bind(c, name = "libsupermesh_integer_map_has_key")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
      integer(kind = c_int) :: present
    end subroutine cinteger_map_has_key

    subroutine cinteger_map_fetch_pair(i, index, key, value) bind(c, name = "libsupermesh_integer_map_fetch_pair")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: index
      integer(kind = c_int) :: key
      integer(kind = c_int) :: value
    end subroutine cinteger_map_fetch_pair
  end interface

  interface allocate
    module procedure allocate_integer_map, allocate_integer_map_rank_1
  end interface allocate
  
  interface deallocate
    module procedure deallocate_integer_map, deallocate_integer_map_rank_1
  end interface deallocate

  interface insert
    module procedure integer_map_insert, integer_map_insert_rank_1
  end interface insert

  interface key_count
    module procedure integer_map_size, integer_map_size_rank_1
  end interface key_count

  interface fetch
    module procedure integer_map_fetch, integer_map_fetch_rank_1
  end interface fetch

  interface remove
    module procedure integer_map_remove, integer_map_remove_rank_1
  end interface remove

  interface has_key
    module procedure integer_map_has_key, integer_map_has_key_rank_1
  end interface has_key

  interface fetch_pair
    module procedure integer_map_fetch_pair, integer_map_fetch_pair_rank_1
  end interface fetch_pair

  public :: integer_map, allocate, deallocate, insert, key_count, &
    & fetch, remove, has_key, fetch_pair

contains 

  subroutine allocate_integer_map(imap)
    type(integer_map), intent(out) :: imap

    allocate(imap%ptr)
    call cinteger_map_new(imap%ptr)

  end subroutine allocate_integer_map
  
  subroutine allocate_integer_map_rank_1(imaps)
    type(integer_map), dimension(:), intent(out) :: imaps
    
    integer :: i
    
    do i = 1, size(imaps)
      allocate(imaps(i)%ptr)
      call cinteger_map_new(imaps(i)%ptr)
    end do
  
  end subroutine allocate_integer_map_rank_1

  subroutine deallocate_integer_map(imap)
    type(integer_map), intent(inout) :: imap

    call cinteger_map_delete(imap%ptr)
    deallocate(imap%ptr)

  end subroutine deallocate_integer_map

  subroutine deallocate_integer_map_rank_1(imaps)
    type(integer_map), dimension(:), intent(inout) :: imaps

    integer :: i

    do i = 1, size(imaps)
      call cinteger_map_delete(imaps(i)%ptr)
      deallocate(imaps(i)%ptr)
    end do

  end subroutine deallocate_integer_map_rank_1

  subroutine integer_map_insert(imap, key, value)
    type(integer_map), intent(inout) :: imap
    integer, intent(in) :: key
    integer, intent(in) :: value

    call cinteger_map_insert(imap%ptr, key, value)

  end subroutine integer_map_insert
  
  subroutine integer_map_insert_rank_1(imap, keys, values)
    type(integer_map), intent(inout) :: imap
    integer, dimension(:), intent(in) :: keys
    integer, dimension(size(keys)), intent(in) :: values
    
    integer :: i
    
    do i = 1, size(keys)
      call cinteger_map_insert(imap%ptr, keys(i), values(i))
    end do
  
  end subroutine integer_map_insert_rank_1      

  function integer_map_size(imap) result(s)
    type(integer_map), intent(inout) :: imap

    integer :: s

    call cinteger_map_size(imap%ptr, s)

  end function integer_map_size
  
  function integer_map_size_rank_1(imaps) result(s)
    type(integer_map), dimension(:), intent(inout) :: imaps
    
    integer, dimension(size(imaps)) :: s
    
    integer :: i
    
    do i = 1, size(imaps)
      call cinteger_map_size(imaps(i)%ptr, s(i))
    end do
    
  end function integer_map_size_rank_1

  function integer_map_fetch(imap, key) result(value)
    type(integer_map), intent(inout) :: imap
    integer, intent(in) :: key

    integer :: value

    call cinteger_map_fetch(imap%ptr, key, value)

  end function integer_map_fetch
  
  function integer_map_fetch_rank_1(imap, keys) result(values)
    type(integer_map), intent(inout) :: imap
    integer, dimension(:), intent(in) :: keys

    integer, dimension(size(keys)) :: values

    integer :: i

    do i = 1, size(keys)
      call cinteger_map_fetch(imap%ptr, keys(i), values(i))
    end do

  end function integer_map_fetch_rank_1

  subroutine integer_map_remove(imap, key)
    type(integer_map), intent(inout) :: imap
    integer, intent(in) :: key

    call cinteger_map_remove(imap%ptr, key)

  end subroutine integer_map_remove
  
  subroutine integer_map_remove_rank_1(imap, keys)
    type(integer_map), intent(inout) :: imap
    integer, dimension(:), intent(in) :: keys
    
    integer :: i
    
    do i = 1, size(keys)
      call cinteger_map_remove(imap%ptr, keys(i))
    end do
  
  end subroutine integer_map_remove_rank_1

  function integer_map_has_key(imap, key) result(present)
    type(integer_map), intent(inout) :: imap
    integer, intent(in) :: key

    logical :: present

    integer(kind = c_int) :: lpresent

    call cinteger_map_has_key(imap%ptr, key, lpresent)
    present = (lpresent /= 0)

  end function integer_map_has_key
  
  function integer_map_has_key_rank_1(imap, keys) result(present)
    type(integer_map), intent(inout) :: imap
    integer, dimension(:), intent(in) :: keys
    
    logical, dimension(size(keys)) :: present
    
    integer :: i
    integer(kind = c_int) :: lpresent
    
    do i = 1, size(keys)
      call cinteger_map_has_key(imap%ptr, keys(i), lpresent)
      present(i) = (lpresent /= 0)
    end do
    
  end function integer_map_has_key_rank_1

  subroutine integer_map_fetch_pair(imap, index, key, value)
    type(integer_map), intent(inout) :: imap
    integer, intent(in) :: index
    integer, intent(out) :: key
    integer, intent(out) :: value

    call cinteger_map_fetch_pair(imap%ptr, index, key, value)

  end subroutine integer_map_fetch_pair

  subroutine integer_map_fetch_pair_rank_1(imap, indices, keys, values)
    type(integer_map), intent(inout) :: imap
    integer, dimension(:), intent(in) :: indices
    integer, dimension(size(indices)), intent(out) :: keys
    integer, dimension(size(indices)), intent(out) :: values

    integer :: i

    do i = 1, size(indices)
      call cinteger_map_fetch_pair(imap%ptr, indices(i), keys(i), values(i))
    end do

  end subroutine integer_map_fetch_pair_rank_1
  

end module libsupermesh_integer_map
