#ifndef NEMOSYS_MESHQUALITYPARAMS_H_
#define NEMOSYS_MESHQUALITYPARAMS_H_

#include "nemosys_export.h"

class NEMOSYS_EXPORT meshQualityParams {
 public:
  virtual ~meshQualityParams() = default;

 protected:
  meshQualityParams() = default;
  meshQualityParams(const meshQualityParams&) = default;
  meshQualityParams(meshQualityParams &&) =default;
  meshQualityParams &operator=(const meshQualityParams &) = default;
  meshQualityParams &operator=(meshQualityParams &&) = default;
};

#endif  // NEMOSYS_MESHQUALITYPARAMS_H_
