#ifndef NEMOSYS_PROTEUSDRIVER_H_
#define NEMOSYS_PROTEUSDRIVER_H_

#include "nemosys_export.h"
#include "NemDriver.H"
#include "proteusHdf5.H"

namespace NEM {
namespace DRV {

/**
 * @brief Driver to parse input JSON file and read Proteus output.
 **/
class NEMOSYS_EXPORT ProteusDriver : public NemDriver {
 public:
  /**
   * @brief Converts Proteus HDF mesh to VTK/Exodus format using input JSON
   * @param fieldFName Proteus format HDF5 field filename
   * @param meshFName Output VTK mesh filename
   * @param edgeSidesetName Name of sideset written to output Exodus file
   * @param exoMeshFName Output Exodus mesh filename
   * @param lowOrder Boolean converting high order cells to low order
   * @param bndryConst Boolean to employ boundary constraint during refinement
   **/
  ProteusDriver(const std::string &_fieldFName, const std::string &_meshFName,
                const std::string &_edgeSidesetName,
                const std::string &_exoMeshFName, bool lowOrder,
                bool bndryConst);

  ~ProteusDriver() override;

  /**
   * @brief high level driver method
   * @param inputjson input JSON file
   * @return ProteusDriver object
   **/
  static ProteusDriver *readJSON(const jsoncons::json &inputjson);

 private:
  std::shared_ptr<proteusHdf5>
      hdf5Obj; /**< hdf5 object used to store Proteus data */
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_PROTEUSDRIVER_H_
