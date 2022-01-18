/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_SMESHUTILS_H_
#define NEMOSYS_SMESHUTILS_H_

#include <SMDS_StdIterator.hxx>

#include <functional>
#include <type_traits>
#include <utility>

namespace NEM {
namespace MSH {

template <class PtrSMDSIterator,
          class VALUE = typename PtrSMDSIterator::element_type::value_type,
          class EqualVALUE = std::equal_to<VALUE>>
struct SM_StdContainerWrapperFromIter {
  PtrSMDSIterator iter;
  using IterT = SMDS_StdIterator<VALUE, PtrSMDSIterator, EqualVALUE>;
  IterT begin() { return {iter}; }
  IterT end() { return {}; }
};

template <class PtrSMDSIterator,
          class VALUE = typename std::decay<
              PtrSMDSIterator>::type::element_type::value_type,
          class EqualVALUE = std::equal_to<VALUE>>
SM_StdContainerWrapperFromIter<typename std::decay<PtrSMDSIterator>::type,
                               VALUE, EqualVALUE>
containerWrapper(PtrSMDSIterator &&iter) {
  return {std::forward<PtrSMDSIterator>(iter)};
}

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_SMESHUTILS_H_
