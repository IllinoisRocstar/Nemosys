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

! The following uses code from
! femtools/Futils.F90 in Fluidity git revision
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

#include "libsupermesh_debug.h"

module libsupermesh_read_triangle

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_precision, only : real_kind

  implicit none

  private

  public :: read_node, read_ele

contains

#if defined __GFORTRAN__ && __GNUC__ == 4 && __GNUC_MINOR__ == 9
  ! Modified version of free_unit in femtools/Futils.F90 in Fluidity git
  ! revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
  function free_unit()
    integer :: free_unit

    logical :: connected

    do free_unit = 100, 1000
       inquire(unit = free_unit, opened = connected)
       if(.not. connected) return
    end do
    libsupermesh_abort("Unable to find I/O unit")

  end function free_unit
#endif

  ! Read Triangle .node file, as described at:
  !   https://www.cs.cmu.edu/~quake/triangle.node.html
  ! Comments and empty lines are not supported, and vertices must be indexed
  ! from one. Trailing lines are ignored.
  !
  ! Expected format differs from that in the Triangle documentation, in that
  ! arbitrary dimensional nodal data is supported.
  subroutine read_node(filename, dim, positions, attributes, boundary_markers)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! dim x nnodes
    real(kind = real_kind), dimension(:, :), allocatable, intent(out) :: &
      & positions
    ! nattrs x nnodes
    real(kind = real_kind), dimension(:, :), allocatable, optional, intent(out) :: &
      & attributes
    ! nbm x nnodes
    integer, dimension(:, :), allocatable, optional, intent(out) :: &
      & boundary_markers

    integer :: i, ind, nnodes, ldim, nattrs, nbm, unit
    integer, dimension(:), allocatable :: boundary_marker
    real(kind = real_kind), dimension(dim) :: coord
    real(kind = real_kind), dimension(:), allocatable :: attribute

#if defined __GFORTRAN__ && __GNUC__ == 4 && __GNUC_MINOR__ == 9
    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")
#else
    open(newunit = unit, file = trim(filename), status = "old", action = "read")
#endif

    read(unit, *) nnodes, ldim, nattrs, nbm
    if(nnodes < 0) then
      libsupermesh_abort("Invalid number of vertices")
    end if
    if(ldim /= dim) then
      libsupermesh_abort("Invalid dimension")
    end if
    if(nattrs < 0) then
      libsupermesh_abort("Invalid number of attributes")
    end if
    if(nbm < 0 .or. nbm > 1) then
      libsupermesh_abort("Invalid number of boundary markers")
    end if

    allocate(positions(dim, nnodes))
    if(present(attributes)) allocate(attributes(nattrs, nnodes))
    if(present(boundary_markers)) allocate(boundary_markers(nbm, nnodes))

    allocate(attribute(nattrs), boundary_marker(nbm))
    do i = 1, nnodes
      if(nattrs > 0) then
        if(nbm > 0) then
          read(unit, *) ind, coord, attribute, boundary_marker
        else
          read(unit, *) ind, coord, attribute
        end if
      else if(nbm > 0) then
        read(unit, *) ind, coord, boundary_marker
      else
        read(unit, *) ind, coord
      end if
      if(i /= ind) then
        libsupermesh_abort("Invalid vertex number")
      end if
      positions(:, i) = coord
      if(present(attributes) .and. nattrs > 0) attributes(:, i) = attribute
      if(present(boundary_markers) .and. nbm > 0) boundary_markers(:, i) = boundary_marker
    end do
    deallocate(attribute, boundary_marker)

    close(unit)

  end subroutine read_node

  ! Read Triangle .ele file, as described at:
  !   https://www.cs.cmu.edu/~quake/triangle.ele.html
  ! Comments and empty lines are not supported, and triangles must be indexed
  ! from one. Trailing lines are ignored.
  !
  ! Expected format differs from that in the Triangle documentation, in that an
  ! arbitrary number of local nodes is permitted.
  subroutine read_ele(filename, dim, enlist, attributes, nnodes)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! loc x nelements
    integer, dimension(:, :), allocatable, intent(out) :: enlist
    ! nattrs x nelements
    real(kind = real_kind), dimension(:, :), allocatable, optional, intent(out) :: &
      & attributes
    integer, optional, intent(in) :: nnodes

    integer :: i, ind, ncell_nodes, nelements, nattrs, unit
    integer, dimension(:), allocatable :: cell_nodes
    real(kind = real_kind), dimension(:), allocatable :: attribute

#if defined __GFORTRAN__ && __GNUC__ == 4 && __GNUC_MINOR__ == 9
    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")
#else
    open(newunit = unit, file = trim(filename), status = "old", action = "read")
#endif

    read(unit, *) nelements, ncell_nodes, nattrs
    if(nelements < 0) then
      libsupermesh_abort("Invalid number of elements")
    end if
    if(nattrs < 0) then
      libsupermesh_abort("Invalid number of attributes")
    end if

    allocate(enlist(ncell_nodes, nelements))
    if(present(attributes)) allocate(attributes(nattrs, nelements))

    allocate(cell_nodes(ncell_nodes), attribute(nattrs))
    do i = 1, nelements
      if(nattrs > 0) then
        read(unit, *) ind, cell_nodes, attribute
      else
        read(unit, *) ind, cell_nodes
      end if
      if(i /= ind) then
        libsupermesh_abort("Invalid element number")
      end if
      if(any(cell_nodes < 1)) then
        libsupermesh_abort("Invalid node")
      end if
      if(present(nnodes)) then
        if(any(cell_nodes > nnodes)) then
          libsupermesh_abort("Invalid node")
        end if
      end if
      enlist(:, i) = cell_nodes
      if(present(attributes) .and. nattrs > 0) attributes(:, i) = attribute
    end do
    deallocate(cell_nodes, attribute)

    close(unit)

  end subroutine read_ele

end module libsupermesh_read_triangle
