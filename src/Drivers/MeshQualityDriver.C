#include <MeshQualityDriver.H>

#include <memory>

#ifdef HAVE_CFMSH
#  include "MeshQuality.H"
#  include "cfmeshQualityParams.H"
#endif

namespace NEM {
namespace DRV {

MeshQualityDriver::MeshQualityDriver(const std::string &_mesh,
                                     const std::string &ofname) {
  // default constructor does standard check mesh process
  // no improvement should be expected
  mesh = meshBase::Create(_mesh);
  mesh->checkMesh(ofname);
  std::cout << "MeshQualityDriver created" << std::endl;
}

MeshQualityDriver::~MeshQualityDriver() {
  delete mesh;
  std::cout << "MeshQualityDriver destroyed" << std::endl;
}

MeshQualityDriver *MeshQualityDriver::readJSON(
    const jsoncons::json &inputjson) {
  MeshQualityDriver *qualdrvobj;
  std::string _mesh = inputjson["Input Mesh File"].as<std::string>();
  std::string ofname = inputjson["Output File"].as<std::string>();
  std::string engine =
      inputjson.get_with_default("Mesh Quality Engine", "default");

  if (!inputjson.contains("Schedule") || engine == "default") {
    // perform a simple check mesh
    qualdrvobj = new MeshQualityDriver(_mesh, ofname);
    return qualdrvobj;
  }

  if (engine == "cfmesh") {
#ifndef HAVE_CFMSH
    std::cerr << "NEMoSys must be recompiled with cfMesh support." << std::endl;
    exit(1);
#else
    auto *params = new cfmshQualityParams();
    auto *mq = new MeshQuality(params);

    // carry out schedules
    for (const auto &jsch : inputjson["Schedule"].array_range()) {
      std::string qImpMtd = jsch["Method"].as<std::string>();
      if (qImpMtd == "meshOptimizer") {
        params->nIterations = jsch["Params"]["NIterations"].as<int>();
        params->nLoops = jsch["Params"]["NLoops"].as<int>();
        params->qualThrsh = jsch["Params"]["QualityThreshold"].as<double>();
        params->nSrfItr = jsch["Params"]["NSurfaceIterations"].as<int>();
        if (jsch["Params"].contains("ConstrainedCellSet")) {
          params->_withConstraint = true;
          params->consCellSet =
              jsch["Params"]["ConstrainedCellSet"].as<std::string>();
        }
        mq->cfmOptimize();
      }
    }

    return nullptr;
#endif
  }

  std::cerr << "Invalid engine specified: " << engine;
  exit(1);
}

}  // namespace DRV
}  // namespace NEM
