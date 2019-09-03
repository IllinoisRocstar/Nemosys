#ifndef NEMOSYS_REFINEDRIVER_H_
#define NEMOSYS_REFINEDRIVER_H_

#include "nemosys_export.h"
#include "NemDriver.H"
#include "meshBase.H"

#include <string>

class NEMOSYS_EXPORT RefineDriver : public NemDriver {
 public:
  RefineDriver(const std::string &_mesh, const std::string &method,
               const std::string &arrayName, double dev_mult, bool maxIsmin,
               double edgescale, const std::string &ofname, bool transferData,
               double sizeFactor = 1.);

  RefineDriver(const std::string &_mesh, const std::string &method,
               double edgescale, const std::string &ofname, bool transferData);

  RefineDriver(const std::string &_mesh, const std::string &method,
               const std::string &arrayName, int order,
               const std::string &ofname, bool transferData);

  static RefineDriver *readJSON(const jsoncons::json &inputjson);
  static RefineDriver *readJSON(const std::string &ifname);

  ~RefineDriver() override;

  meshBase *getMesh() const { return mesh; }

 private:
  meshBase *mesh;
};

#endif  // NEMOSYS_REFINEDRIVER_H_
