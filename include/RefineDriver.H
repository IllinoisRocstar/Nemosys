#ifndef NEMOSYS_REFINEDRIVER_H_
#define NEMOSYS_REFINEDRIVER_H_

#include <string>

#include "nemosys_export.h"
#include "NemDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT RefineDriver : public NemDriver {
 public:
  RefineDriver() = default;

  RefineDriver(const std::string &_mesh, const std::string &method,
               const std::string &arrayName, double dev_mult, bool maxIsmin,
               double edgescale, const std::string &ofname, bool transferData,
               double sizeFactor = 1.);

  RefineDriver(const std::string &_mesh, const std::string &method,
               double edgescale, const std::string &ofname, bool transferData);

  RefineDriver(const std::string &_mesh, const std::string &method,
               const std::string &arrayName, int order,
               const std::string &ofname, bool transferData);

#ifdef HAVE_CFMSH
  /** @brief This method lets users input a size field or simulation quantity
              field defined over mesh and perform refinement/derefinement.
              It supports one step refinement/unrefinement as well as online
              refinement parallel to solver simulation
      @param _mesh Input mesh file name with extension. If using OpenFOAM mesh,
             define this as empty string (i.e "")
      @param ofname Output mesh file name
      @param method Refinement method (here, "AMR")
      @param inputFile Input field for refinement.
      @param refineInterval Refinement interval/Frequency
      @param lowerRefineLevel Lower refinement level based on field
      @param upperRefineLevel Upper refinement level based on field
      @param unrefineAbove Criteria for unrefinement above this value
      @param unrefineBelow Criteria for unrefinement below this value
      @param writeFieldData A boolean to enable writing of updated field data
      @param writeMesh A boolean to enable writing mesh at each refine interval
      @param writeRefHistory A boolean to enable writing of refinement history
      @param timeStep If running with simulation. Default is 1.
      @param endTime If running with simulation. Default is 1.
      @param nBufferLayers If user wants slower than 2:1 refinement.
      @param maxCells Maximum allowed cells in mesh. Once this limit is reached,
             refinement stops.
      @param refCriteria Refinement based on "Value" or "Gradient".
      @param startT Simulation start time.
  **/
  RefineDriver(const std::string &_mesh, const std::string &ofname,
               const std::string &method, const std::string &inputFile,
               const int &refineInterval, const int &maxRefinement,
               const double &lowerRefineLevel, const double &upperRefineLevel,
               const double &unrefineAbove, const double &unrefineBelow,
               const bool &writeFieldData, const bool &writeMesh,
               const bool &writeRefHistory, const double &timeStep,
               const double &endTime, const int &nBufferLayers, int &maxCells,
               const std::string &refCriteria, const double &startT,
               const std::string &MLName);
#endif

  static RefineDriver *readJSON(const jsoncons::json &inputjson);
  static RefineDriver *readJSON(const std::string &ifname);

  ~RefineDriver() override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_REFINEDRIVER_H_
