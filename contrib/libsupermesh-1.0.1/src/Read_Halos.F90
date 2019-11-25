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

! The following code is derived from femtools/Halo_Data_Types.F90 and
! femtools/Halos_Registration.F90 in Fluidity git revision
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

module libsupermesh_read_halos

  use iso_c_binding, only : c_char, c_int, c_null_char, c_ptr
  use iso_fortran_env, only : real64
  use mpi

  use libsupermesh_debug, only : abort_pinpoint

  implicit none

  private

  public :: int_array, halo_type, read_halo, deallocate

  interface
    subroutine cread_halo(data, filename, process, nprocs) bind(c, name = "libsupermesh_read_halo")
      use iso_c_binding, only : c_char, c_int, c_ptr
      implicit none
      type(c_ptr) :: data
      character(kind = c_char), dimension(*) :: filename
      integer(kind = c_int) :: process
      integer(kind = c_int) :: nprocs
    end subroutine cread_halo

    subroutine chalo_sizes(data, level, nsends, nreceives) bind(c, name = "libsupermesh_halo_sizes")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: data
      integer(kind = c_int), value :: level
      integer(kind = c_int), dimension(*) :: nsends
      integer(kind = c_int), dimension(*) :: nreceives
    end subroutine chalo_sizes

    subroutine chalo_data(data, level, npnodes, send, recv) bind(c, name = "libsupermesh_halo_data")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: data
      integer(kind = c_int), value :: level
      integer(kind = c_int) :: npnodes
      integer(kind = c_int), dimension(*) :: send      
      integer(kind = c_int), dimension(*) :: recv
    end subroutine chalo_data
    
    subroutine cdeallocate_halo(data) bind(c, name = "libsupermesh_deallocate_halo")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: data
    end subroutine cdeallocate_halo
  end interface

  type int_array
    integer, dimension(:), pointer :: v
  end type int_array

  type halo_type
    integer :: comm
    integer :: level
    integer :: process
    integer :: nprocs
    integer :: npnodes
    type(int_array), dimension(:), pointer :: send
    type(int_array), dimension(:), pointer :: recv
  end type halo_type

  interface deallocate
    module procedure deallocate_halo
  end interface deallocate

contains

  subroutine read_halo(basename, halo, level, comm)
    character(len = *), intent(in) :: basename
    type(halo_type), intent(out) :: halo
    integer, optional, intent(in) :: level
    integer, optional, intent(in) :: comm

    character(len = len_trim(basename) + int(log10(real(huge(halo%process), kind = real64))) + 8) :: &
      & filename
    integer :: i, ierr, index
    integer(kind = c_int) :: nprocs, process
    integer(kind = c_int), dimension(:), allocatable :: nsends, nreceives, &
      & send, recv
    type(c_ptr) :: data
    
    if(present(comm)) then
      halo%comm = comm
    else
      halo%comm = MPI_COMM_WORLD
    end if

    call MPI_Comm_rank(halo%comm, halo%process, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Comm_size(halo%comm, halo%nprocs, ierr);  assert(ierr == MPI_SUCCESS)        
    write(filename, "(a,a,i0,a)") trim(basename), "_", halo%process, ".halo"
    call cread_halo(data, to_cstring(filename), process, nprocs)
    if(nprocs /= halo%nprocs) then
      libsupermesh_abort("Failed to read halo file '" // trim(filename) // "': Unexpected number of processes")
    else if(process /= halo%process) then
      libsupermesh_abort("Failed to read halo file '" // trim(filename) // "': Unexpected process number")
    end if

    if(present(level)) then
      halo%level = level
    else
      halo%level = 2
    end if
    allocate(nsends(halo%nprocs), nreceives(halo%nprocs))
    call chalo_sizes(data, halo%level, nsends, nreceives)
    
    allocate(halo%send(halo%nprocs), halo%recv(halo%nprocs))
    do i = 1, halo%nprocs
      allocate(halo%send(i)%v(nsends(i)), halo%recv(i)%v(nreceives(i)))
    end do
    allocate(send(sum(nsends)), recv(sum(nreceives)))
    call chalo_data(data, halo%level, halo%npnodes, send, recv)
    index = 1
    do i = 1, halo%nprocs
      halo%send(i)%v = send(index:index + nsends(i) - 1)
      index = index + nsends(i)
    end do
    index = 1
    do i = 1, halo%nprocs
      halo%recv(i)%v = recv(index:index + nreceives(i) - 1)
      index = index + nreceives(i)
    end do    
    deallocate(nsends, nreceives, send, recv)
    
    call cdeallocate_halo(data)

  end subroutine read_halo
    
  pure function to_cstring(fs) result(cs)
    character(len = *), intent(in) :: fs
    
    character(kind = c_char), dimension(len_trim(fs) + 1) :: cs
    
    integer :: i
    
    forall(i = 1:size(cs) - 1)
      cs(i) = fs(i:i)
    end forall
    cs(size(cs)) = c_null_char
    
  end function to_cstring

  pure elemental subroutine deallocate_halo(halo)
    type(halo_type), intent(inout) :: halo

    integer :: i

    do i = 1, halo%nprocs
      deallocate(halo%send(i)%v, halo%recv(i)%v)
    end do
    deallocate(halo%send, halo%recv)

  end subroutine deallocate_halo

end module libsupermesh_read_halos
