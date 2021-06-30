#ifndef MESHINGPARAMS_H
#define MESHINGPARAMS_H

#include "nemosys_export.h"

class NEMOSYS_EXPORT meshingParams {
 public:
  virtual ~meshingParams() = default;

 protected:
  meshingParams() = default;
  meshingParams(const meshingParams &) = default;
  meshingParams(meshingParams &&) = default;
  meshingParams &operator=(const meshingParams &) = default;
  meshingParams &operator=(meshingParams &&) = default;
};

#endif
