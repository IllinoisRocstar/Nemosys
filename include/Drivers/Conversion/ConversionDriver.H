#ifndef NEMOSYS_CONVERSIONDRIVER_H_
#define NEMOSYS_CONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

#include <map>
#include <string>
#include <vector>

#include "Mesh/exoMesh.H"
#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT ConversionDriver : public NemDriver {
 public:
  static void genExo(std::vector<meshBase *> meshes, const std::string &fname);

  JSONCONS_TYPE_TRAITS_FRIEND

 protected:
  static void genExo(meshBase *mb, NEM::MSH::EXOMesh::exoMesh *em,
                     const int &ndeIdOffset, const int &elmIdOffset, int &ins,
                     int &ieb, int &iss, std::string mshName,
                     const bool &usePhys, int &ndeIdOffset_local,
                     int &elmIdOffset_local, const bool &makeFreeSurfSS,
                     const bool &splitTopBotSS,
                     std::vector<std::string> sideSetNames);

  static void procExo(const jsoncons::json &ppJson, const std::string &fname,
                      NEM::MSH::EXOMesh::exoMesh *em);

  /**
   * @brief Creates side set(s) for the free surface, exterior surface, during
   * conversion
   * @param mb Meshbase object
   * @param em ExoMesh object
   * @param elmIdOffset Exodus element ID offset
   * @param v2e_elemID_map VTK to EXO element ID map
   * @param splitTopBotSS Boolean to split the side set into three side sets
   * @param sideSetNames The list of name(s) for the side set(s)
   */
  static void freeSurfaceSideSet(const meshBase *mb,
                                 NEM::MSH::EXOMesh::exoMesh *em,
                                 int elmIdOffset,
                                 std::map<int, int> v2e_elemID_map,
                                 bool splitTopBotSS,
                                 std::vector<std::string> sideSetNames);

 private:
  static constexpr const char *programType = "Conversion";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_CONVERSIONDRIVER_H_
