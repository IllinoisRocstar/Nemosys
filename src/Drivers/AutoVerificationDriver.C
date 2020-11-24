#include "AutoVerificationDriver.H"

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace DRV {
AutoVerificationDriver::AutoVerificationDriver(
    meshBase *coarseMesh, meshBase *fineMesh, meshBase *finerMesh,
    std::vector<int> arrayIDs, std::string transferType, double targetGCI,
    int numThreads) {
  std::cout << "AutoVerificationDriver created." << std::endl;
#ifdef HAVE_OPENMP
  omp_set_num_threads(numThreads);
  std::cout << "Number of threads set to : " << numThreads << std::endl;
#endif
  std::cout << "Running verification." << std::endl;
  this->oac = std::make_shared<OrderOfAccuracy>(OrderOfAccuracy(
      coarseMesh, fineMesh, finerMesh, arrayIDs, transferType, targetGCI));
  std::cout << "Checking if in asymptotic range." << std::endl;
  std::cout << "Target GCI is set to : " << oac->getTargetGCI() << std::endl;
  auto asymp = this->oac->checkAsymptoticRange();
  bool inRange = true;
  for (int i = 0; i < asymp.size(); ++i) {
    for (int j = 0; j < asymp[0].size(); ++j) {
      double gci = asymp[i][j];
      if (gci > oac->getTargetGCI()) {
        std::cout << "GCI of " << gci;
        std::cout << " exceeds target GCI of " << targetGCI << std::endl;
        std::cout << "at array " << i << ", component " << j << std::endl;
        inRange = false;
      }
    }
  }
  if (inRange) {
    std::cout << "Grid is in target asymptotic range." << std::endl;
  } else {
    std::cout << "Grid is not in target asymptotic range." << std::endl;
  }
}

AutoVerificationDriver *AutoVerificationDriver::readJSON(
    const std::string &ifname) {
  if (nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  std::ifstream inputStream(ifname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array()) {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed\n";
    return AutoVerificationDriver::readJSON(inputjson[0]);
  } else {
    return AutoVerificationDriver::readJSON(inputjson);
  }
}

AutoVerificationDriver *AutoVerificationDriver::readJSON(
    const jsoncons::json &inputjson) {
  auto meshFileOptions = inputjson["Mesh File Options"];

  std::string coarseMeshFname =
      meshFileOptions["Coarse Mesh File"].as<std::string>();
  std::string fineMeshFname =
      meshFileOptions["Fine Mesh File"].as<std::string>();
  std::string finerMeshFname =
      meshFileOptions["Finer Mesh File"].as<std::string>();

  auto coarseMesh = meshBase::CreateShared(coarseMeshFname);
  auto fineMesh = meshBase::CreateShared(fineMeshFname);
  auto finerMesh = meshBase::CreateShared(finerMeshFname);

  auto verificationOptions = inputjson["Verification Options"];

  // required
  std::vector<int> arrayIds;
  jsoncons::json jsonArrayIDs = inputjson["Verification Options"]["Array IDs"];
  if (!jsonArrayIDs.is_null()) {
    for (const auto &item : jsonArrayIDs.array_range()) {
      int id = item.as<int>();
      arrayIds.push_back(id);
    }
  } else {
    std::cerr << "No array IDs specified for verification. Aborting."
              << std::endl;
    exit(1);
  }

  // optional
  std::string transferType = "Consistent Interpolation";
  if (verificationOptions.find("Transfer Type") !=
      verificationOptions.object_range().end()) {
    transferType = verificationOptions["Transfer Type"].as<std::string>();
  }
  // optional
  double targetGCI = 1.1;
  if (verificationOptions.find("Target GCI") !=
      verificationOptions.object_range().end()) {
    targetGCI = verificationOptions["Target GCI"].as<double>();
  }

  int numThreads = 1;
#ifdef HAVE_OPENMP
  numThreads = omp_get_max_threads();
  if (verificationOptions.find("Threads") !=
      verificationOptions.object_range().end()) {
    numThreads = verificationOptions["Threads"].as<int>();
  }
#else
  if (verificationOptions.find("Threads") !=
      verificationOptions.object_range().end()) {
    std::cerr << "OpenMP is not enabled. Verification will continue in serial."
              << std::endl;
  }
#endif
  return new AutoVerificationDriver(coarseMesh.get(), fineMesh.get(),
                                    finerMesh.get(), arrayIds, transferType,
                                    targetGCI, numThreads);
}
}  // namespace DRV
}  // namespace NEM
