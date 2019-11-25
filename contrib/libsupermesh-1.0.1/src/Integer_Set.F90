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
! femtools/Integer_set.F90 in Fluidity git revision
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

module libsupermesh_integer_set

  use iso_c_binding, only : c_int, c_ptr
  
  implicit none
  
  private
  
  type integer_set
    type(c_ptr), pointer :: ptr
  end type integer_set

  interface
    subroutine cinteger_set_new(i) bind(c, name = "libsupermesh_integer_set_new")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine cinteger_set_new

    subroutine cinteger_set_delete(i) bind(c, name = "libsupermesh_integer_set_delete")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine cinteger_set_delete

    subroutine cinteger_set_insert(i, value, changed) bind(c, name = "libsupermesh_integer_set_insert")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: value
      integer(kind = c_int) :: changed
    end subroutine cinteger_set_insert

    subroutine cinteger_set_size(i, size) bind(c, name = "libsupermesh_integer_set_size")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: size
    end subroutine cinteger_set_size

    subroutine cinteger_set_fetch(i, index, value) bind(c, name = "libsupermesh_integer_set_fetch")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: index
      integer(kind = c_int) :: value
    end subroutine cinteger_set_fetch

    subroutine cinteger_set_remove(i, value) bind(c, name = "libsupermesh_integer_set_remove")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: value
    end subroutine cinteger_set_remove

    subroutine cinteger_set_has_value(i, value, present) bind(c, name = "libsupermesh_integer_set_has_value")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: value
      integer(kind = c_int) :: present
    end subroutine cinteger_set_has_value
  end interface

  interface allocate
    module procedure allocate_integer_set, allocate_integer_set_rank_1
  end interface allocate
  
  interface deallocate
    module procedure deallocate_integer_set, deallocate_integer_set_rank_1
  end interface deallocate

  interface insert
    module procedure integer_set_insert, integer_set_insert_rank_1
  end interface insert
  
  interface key_count
    module procedure integer_set_size, integer_set_size_rank_1
  end interface key_count
  
  interface fetch
    module procedure integer_set_fetch, integer_set_fetch_rank_1
  end interface fetch

  interface remove
    module procedure integer_set_remove, integer_set_remove_rank_1
  end interface remove

  interface has_value
    module procedure integer_set_has_value, integer_set_has_value_rank_1
  end interface has_value

  public :: integer_set, allocate, deallocate, insert, key_count, fetch, &
    & remove, has_value

contains 
  
  subroutine allocate_integer_set(iset)
    type(integer_set), intent(out) :: iset

    allocate(iset%ptr)
    call cinteger_set_new(iset%ptr)

  end subroutine allocate_integer_set
  
  subroutine allocate_integer_set_rank_1(isets)
    type(integer_set), dimension(:), intent(out) :: isets
    
    integer :: i
    
    do i = 1, size(isets)
      allocate(isets(i)%ptr)
      call cinteger_set_new(isets(i)%ptr)
    end do
  
  end subroutine allocate_integer_set_rank_1

  subroutine deallocate_integer_set(iset)
    type(integer_set), intent(inout) :: iset

    call cinteger_set_delete(iset%ptr)
    deallocate(iset%ptr)

  end subroutine deallocate_integer_set

  subroutine deallocate_integer_set_rank_1(isets)
    type(integer_set), dimension(:), intent(inout) :: isets
    
    integer :: i
    
    do i = 1, size(isets)
      call cinteger_set_delete(isets(i)%ptr)
      deallocate(isets(i)%ptr)
    end do
    
  end subroutine deallocate_integer_set_rank_1

  subroutine integer_set_insert(iset, value, changed)
    type(integer_set), intent(inout) :: iset
    integer, intent(in) :: value
    logical, intent(out), optional :: changed

    integer(kind = c_int) :: lchanged

    call cinteger_set_insert(iset%ptr, value, lchanged)
    if(present(changed)) changed = (lchanged /= 0)

  end subroutine integer_set_insert

  subroutine integer_set_insert_rank_1(iset, values)
    type(integer_set), intent(inout) :: iset
    integer, dimension(:), intent(in) :: values

    integer :: i
    integer(kind = c_int) :: lchanged

    do i = 1, size(values)
      call cinteger_set_insert(iset%ptr, values(i), lchanged)
    end do

  end subroutine integer_set_insert_rank_1
  
  function integer_set_size(iset) result(s)
    type(integer_set), intent(inout) :: iset

    integer :: s

    call cinteger_set_size(iset%ptr, s)

  end function integer_set_size
  
  function integer_set_size_rank_1(isets) result(s)
    type(integer_set), dimension(:), intent(inout) :: isets

    integer, dimension(size(isets)) :: s
    
    integer :: i
    
    do i = 1, size(isets)
      call cinteger_set_size(isets(i)%ptr, s(i))
    end do
  
  end function integer_set_size_rank_1

  function integer_set_fetch(iset, index) result(value)
    type(integer_set), intent(inout) :: iset
    integer, intent(in) :: index

    integer :: value

    call cinteger_set_fetch(iset%ptr, index, value)

  end function integer_set_fetch

  function integer_set_fetch_rank_1(iset, indices) result(values)
    type(integer_set), intent(inout) :: iset
    integer, dimension(:), intent(in) :: indices

    integer, dimension(size(indices)) :: values
    
    integer :: i
    
    do i = 1, size(indices)
      call cinteger_set_fetch(iset%ptr, indices(i), values(i))
    end do

  end function integer_set_fetch_rank_1

  subroutine integer_set_remove(iset, value)
    type(integer_set), intent(inout) :: iset
    integer, intent(in) :: value

    call cinteger_set_remove(iset%ptr, value)
  
  end subroutine integer_set_remove

  subroutine integer_set_remove_rank_1(iset, values)
    type(integer_set), intent(inout) :: iset
    integer, dimension(:), intent(in) :: values

    integer :: i

    do i = 1, size(values)
      call cinteger_set_remove(iset%ptr, values(i))
    end do
  
  end subroutine integer_set_remove_rank_1

  function integer_set_has_value(iset, value) result(present)
    type(integer_set), intent(inout) :: iset
    integer, intent(in) :: value

    logical :: present

    integer(kind = c_int) :: lpresent
    
    call cinteger_set_has_value(iset%ptr, value, lpresent)
    present = (lpresent /= 0)

  end function integer_set_has_value

  function integer_set_has_value_rank_1(iset, values) result(present)
    type(integer_set), intent(inout) :: iset
    integer, dimension(:), intent(in) :: values

    logical, dimension(size(values)) :: present
    
    integer:: i
    integer(kind = c_int) :: lpresent
    
    do i = 1, size(values)
      call cinteger_set_has_value(iset%ptr, values(i), lpresent)
      present(i) = (lpresent /= 0)
    end do

  end function integer_set_has_value_rank_1

end module libsupermesh_integer_set
