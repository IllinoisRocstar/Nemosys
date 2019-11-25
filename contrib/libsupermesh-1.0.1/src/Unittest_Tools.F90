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
! femtools/Unittest_tools.F90 in Fluidity git revision
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

module libsupermesh_unittest_tools
  ! Utility functions for unit testing

  use iso_fortran_env, only : output_unit
  
  use libsupermesh_precision, only : real_kind

  implicit none

  private
  
  public :: report_test, operator(.feq.), operator(.fne.), fequals, fnequals

  interface operator(.feq.)
    module procedure fequals_op, fequals_array_scalar_op, &
      & fequals_matrix_scalar_op
  end interface operator(.feq.)

  interface operator(.fne.)
    module procedure fnequals_op, fnequals_array_scalar_op, &
      & fnequals_matrix_scalar_op
  end interface operator(.fne.)

  interface fequals
    module procedure fequals, fequals_array_scalar, fequals_matrix_scalar
  end interface fequals

  interface fnequals
    module procedure fnequals, fnequals_array_scalar, fnequals_matrix_scalar 
  end interface fnequals

contains

  subroutine report_test(title, fail, warn, msg, unit)
    !!< This is the subroutine used by unit tests to report the output of
    !!< a test case.

    !! title: the name of the test case.
    character(len = *), intent(in) :: title
    !! msg: an explanatory message printed if the test case fails.
    character(len = *), intent(in) :: msg
    !! Has the test case failed, or triggered a warning? Set fail or warn to .true. if so.
    logical, intent(in) :: fail, warn
    integer, optional, intent(in) :: unit
    
    integer :: lunit
    
    if(present(unit)) then
      lunit = unit
    else
      lunit = output_unit
    end if

    if(fail) then
      write(lunit, "(a)") "(Fail: " // trim(title) // "; error: " // trim(msg) // ")"
    else if(warn) then
      write(lunit, "(a)") "(Warn: " // trim(title) // "; error: " // trim(msg) // ")"
    else
      write(lunit, "(a)") "(Pass: " // trim(title) // ")"
    end if

  end subroutine report_test

  pure elemental function fequals_op(float1, float2) result(equals)
    real(kind = real_kind), intent(in) :: float1
    real(kind = real_kind), intent(in) :: float2

    logical :: equals

    equals = fequals(float1, float2)

  end function fequals_op

  pure elemental function fnequals_op(float1, float2) result(nequals)
    real(kind = real_kind), intent(in) :: float1
    real(kind = real_kind), intent(in) :: float2

    logical :: nequals

    nequals = fnequals(float1, float2)

  end function fnequals_op

  pure function fequals_array_scalar_op(array1, float2) result(equals)
    real(kind = real_kind), dimension(:), intent(in) :: array1
    real(kind = real_kind), intent(in) :: float2

    logical :: equals

    equals = fequals(array1, float2)

  end function fequals_array_scalar_op

  pure function fnequals_array_scalar_op(array1, float2) result(nequals)
    real(kind = real_kind), dimension(:), intent(in) :: array1
    real(kind = real_kind), intent(in) :: float2

    logical :: nequals

    nequals = fnequals(array1, float2)

  end function fnequals_array_scalar_op
  
  pure function fequals_matrix_scalar_op(mat1, float2) result(equals)
    real(kind = real_kind), dimension(:, :), intent(in) :: mat1
    real(kind = real_kind), intent(in) :: float2
    
    logical :: equals

    equals = fequals(mat1, float2)

  end function fequals_matrix_scalar_op
  
  pure function fnequals_matrix_scalar_op(mat1, float2) result(nequals)
    real(kind = real_kind), dimension(:, :), intent(in) :: mat1
    real(kind = real_kind), intent(in) :: float2
    
    logical :: nequals

    nequals = fnequals(mat1, float2)

  end function fnequals_matrix_scalar_op

  pure elemental function fequals(float1, float2, tol) result(equals)
    real(kind = real_kind), intent(in) :: float1
    real(kind = real_kind), intent(in) :: float2
    real(kind = real_kind), intent(in), optional :: tol
    
    logical :: equals

    real(kind = real_kind) :: eps
    
    if(present(tol)) then
      eps = abs(tol)
    else
      eps = 1.0e2_real_kind * max(epsilon(eps), spacing(float1), spacing(float2))
    end if
    equals = abs(float1 - float2) < eps

  end function fequals

  pure elemental function fnequals(float1, float2, tol) result(nequals)
    real(kind = real_kind), intent(in) :: float1
    real(kind = real_kind), intent(in) :: float2
    real(kind = real_kind), optional, intent(in) :: tol

    logical :: nequals

    nequals = .not. fequals(float1, float2, tol = tol)

  end function fnequals

  pure function fequals_array_scalar(array1, float2, tol) result(equals)
    real(kind = real_kind), dimension(:), intent(in) :: array1
    real(kind = real_kind), intent(in) :: float2
    real(kind = real_kind), intent(in), optional :: tol

    logical :: equals
    
    integer :: i

    do i = 1, size(array1)
      if(fnequals(array1(i), float2, tol = tol)) then
        equals = .false.
        return
      end if
    end do
    equals = .true.

  end function fequals_array_scalar

  pure function fnequals_array_scalar(array1, float2, tol) result(nequals)
    real(kind = real_kind), dimension(:), intent(in) :: array1
    real(kind = real_kind), intent(in) :: float2
    real(kind = real_kind), intent(in), optional :: tol

    logical :: nequals

    nequals = .not. fequals(array1, float2, tol = tol)

  end function fnequals_array_scalar

  pure function fequals_matrix_scalar(mat1, float2, tol) result(equals)
    real(kind = real_kind), dimension(:, :), intent(in) :: mat1
    real(kind = real_kind), intent(in) :: float2
    real(kind = real_kind), optional, intent(in) :: tol
    
    logical :: equals

    integer :: i

    do i = 1, size(mat1, 1)
      if(fnequals(mat1(i, :), float2, tol = tol)) then
        equals = .false.
        return
      end if
    end do
    equals = .true.

  end function fequals_matrix_scalar
  
  pure function fnequals_matrix_scalar(mat1, float2, tol) result(nequals)
    real(kind = real_kind), dimension(:, :), intent(in) :: mat1
    real(kind = real_kind), intent(in) :: float2
    real(kind = real_kind), optional, intent(in) :: tol
    
    logical :: nequals

    nequals = .not. fequals(mat1, float2, tol = tol)

  end function fnequals_matrix_scalar

end module libsupermesh_unittest_tools
