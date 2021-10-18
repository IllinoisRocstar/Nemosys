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
