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

#include "libsupermesh_debug.h"

module libsupermesh

  use libsupermesh_intersection_finder
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_precision, only : real_kind
  use libsupermesh_supermesh

  implicit none

  public
  
  integer, parameter :: version_major = LIBSUPERMESH_VERSION_MAJOR
  integer, parameter :: version_minor = LIBSUPERMESH_VERSION_MINOR
  integer, parameter :: version_patch = LIBSUPERMESH_VERSION_PATCH
#ifdef LIBSUPERMESH_VERSION_RELEASE
  logical, parameter :: version_release = .true.
#else
  logical, parameter :: version_release = .false.
#endif
  character(len = *), parameter :: version = LIBSUPERMESH_VERSION

end module libsupermesh
