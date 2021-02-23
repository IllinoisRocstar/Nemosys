#ifndef NEMOSYS_CONVERSIONDRIVER_H_
#define NEMOSYS_CONVERSIONDRIVER_H_

#include <string>

#include "nemosys_export.h"
#include "NemDriver.H"
#include "exoMesh.H"
#include "meshBase.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT ConversionDriver : public NemDriver {
 public:
  ConversionDriver() : source(nullptr) {}

  /**
   * @brief ConversionDriver constructor
   * @param srcmsh Source mesh file name
   * @param trgmsh Target mesh file name
   * @param method Conversion method
   * @param ofname Output file name
   * @param inputjson Input json file
   */
  ConversionDriver(const std::string &srcmsh, const std::string &trgmsh,
                   const std::string &method, const std::string &ofname,
                   const jsoncons::json &inputjson);
  ~ConversionDriver() override;

  /**
   * @brief Reads the input json file
   * @param inputjson Input json object
   * @return ConversionDriver object
   */
  static ConversionDriver *readJSON(const jsoncons::json &inputjson);

  /**
   * @brief Reads the input json file
   * @param inputjson Input json object
   * @return ConversionDriver object
   */
  static ConversionDriver *readJSON(const std::string &ifname);

 public:
  /**
   * @brief Generates an EXODUS database
   * @param meshes Vector of meshBase objects (meshes)
   * @param fname Output file name
   */
  static void genExo(std::vector<meshBase *> meshes, const std::string &fname);

 private:
  /**
   * @brief Generates an EXODUS database
   * @param opts The 'Conversion Options' json object
   * @param fname Output file name
   */
  static void genExo(const jsoncons::json &opts, const std::string &fname);

  /**
   * @brief Generates an EXODUS database
   * @param mb Meshbase object
   * @param em ExoMesh object
   * @param ndeIdOffset Exodus node ID offset
   * @param elmIdOffset Exodus element ID offset
   * @param ins nodeset counter
   * @param ieb element block counter
   * @param iss sideset counter
   * @param mshName Prefix name for element blocks and global nodeset
   * @param usePhys Boolean to use Gmsh physical groups for element blocks
   * @param ndeIdOffset_local Exodus node ID offset per input msh file
   * @param elmIdOffset_local Exodus element block ID offset per input msh file
   * @param makeFreeSurfSS Make sideset(s) on free surfaces (exterior surfaces)
   * @param splitTopBotSS Boolean to split the side set into three side sets
   * @param sideSetNames The list of name(s) for the side set(s)
   */
  static void genExo(meshBase *mb, NEM::MSH::EXOMesh::exoMesh *em,
                     const int &ndeIdOffset, const int &elmIdOffset, int &ins,
                     int &ieb, int &iss, std::string mshName,
                     const bool &usePhys, int &ndeIdOffset_local,
                     int &elmIdOffset_local, const bool &makeFreeSurfSS,
                     const bool &splitTopBotSS,
                     std::vector<std::string> sideSetNames);
  /**
   * @brief Performs post processing on EXODUS database
   * @param ppJson Post Processing section of input json
   * @param fname Output EXODUS mesh file name
   * @param em EXODUS database object
   */
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

  /**
   * @brief Manipulates an existing EXODUSII mesh
   * @param opts The 'Conversion Options' json object
   * @param fname Output file name
   */
  static void manipExo(const jsoncons::json &opts, const std::string &fname);

 private:
  meshBase *source;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_CONVERSIONDRIVER_H_
