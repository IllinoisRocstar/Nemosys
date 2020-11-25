#include <ProteusDriver.H>
#include <meshBase.H>
#include <proteusHdf5.H>

namespace NEM {
namespace DRV {

ProteusDriver::ProteusDriver(const std::string &fieldFName,
                             const std::string &meshFName,
                             const std::string &edgeSidesetName,
                             const std::string &exoMeshFName, bool lowOrder,
                             bool bndryConst) {
  std::cout << "ProteusDriver created\n";

  // Create HDF5 class of Proteus type
  hdf5Obj =
      std::make_shared<proteusHdf5>(fieldFName, meshFName, edgeSidesetName,
                                    exoMeshFName, lowOrder, bndryConst);
}

// Parse JSON input file
ProteusDriver *ProteusDriver::readJSON(const jsoncons::json &inputjson) {
  std::string fieldFName = inputjson["HDF5 Field File"].as<std::string>();
  std::string meshFName = inputjson["Output Mesh File"].as<std::string>();
  bool lowOrder = inputjson.get_with_default("Output Low Order", false);
  std::string edgeSidesetName = inputjson["Edge Sideset"].as<std::string>();
  std::string exoMeshFName = inputjson["Output Exodus File"].as<std::string>();
  bool bndryConst = inputjson.get_with_default("Constraint Boundary", true);

  return new ProteusDriver(fieldFName, meshFName, edgeSidesetName, exoMeshFName,
                           lowOrder, bndryConst);
}

// Destructor
ProteusDriver::~ProteusDriver() {
  std::cout << "ProteusDriver destroyed" << std::endl;
}

}  // namespace DRV
}  // namespace NEM
