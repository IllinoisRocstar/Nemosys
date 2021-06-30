#include "Drivers/Conversion/VtkToPatranConversionDriver.H"

#include "patran.H"

namespace NEM {
namespace DRV {

VtkToPatranConversionDriver::BoundaryCond::BoundaryCond(int patchNum)
    : patchNum(patchNum) {}

VtkToPatranConversionDriver::FaceBC::FaceBC(int patchNum, int rocFracFSIType)
    : BoundaryCond(patchNum), rocFracFSIType(rocFracFSIType) {}

jsoncons::string_view VtkToPatranConversionDriver::FaceBC::getBCType() const {
  return bcType;
}

VtkToPatranConversionDriver::NodeBC::NodeBC(int patchNum,
                                            int rocfracControlType,
                                            bool structural, bool meshMotion,
                                            bool thermal)
    : BoundaryCond(patchNum),
      rocfracControlType(rocfracControlType),
      structural(structural),
      meshMotion(meshMotion),
      thermal(thermal) {}

jsoncons::string_view VtkToPatranConversionDriver::NodeBC::getBCType() const {
  return bcType;
}

VtkToPatranConversionDriver::Opts::Opts(
    std::vector<std::shared_ptr<BoundaryCond>> bcInfo,
    std::vector<int> nodePatchPreference)
    : bcInfo(std::move(bcInfo)),
      nodePatchPreference(std::move(nodePatchPreference)) {}

VtkToPatranConversionDriver::VtkToPatranConversionDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

VtkToPatranConversionDriver::VtkToPatranConversionDriver()
    : VtkToPatranConversionDriver({{}, {}}, {{}, {}}) {}

const VtkToPatranConversionDriver::Files &
VtkToPatranConversionDriver::getFiles() const {
  return files_;
}

void VtkToPatranConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const VtkToPatranConversionDriver::Opts &VtkToPatranConversionDriver::getOpts()
    const {
  return opts_;
}

void VtkToPatranConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void VtkToPatranConversionDriver::execute() const {
  if (this->files_.inputMeshFile.find(".vt") != std::string::npos) {
    std::cout << "Detected file in VTK format" << std::endl;
    std::cout << "Converting to PATRAN ...." << std::endl;
  } else {
    std::cerr << "Source mesh file is not in VTK format" << std::endl;
  }

  std::ifstream meshStream(this->files_.inputMeshFile);
  if (!meshStream.good()) {
    std::cerr << "Error opening file " << this->files_.inputMeshFile
              << std::endl;
    exit(1);
  }

  // looping through blocks
  std::cout << "Number of Boundary Conditions read: "
            << this->opts_.bcInfo.size() << std::endl;

  // PATRAN specific BC information
  // Each map below maps the indicated information from the "Patch Number"
  // specified in the file

  std::map<int, int> faceTypeMap;  // Rocfrac FSI Type;
  // 0 = no FSI, 1 = FSI w/ burn, 2 = FSI w/o burn, etc.
  // see Rocfrac manual for details
  std::map<int, int>
      nodeTypeMap;  // patch numbers as specified in RocfracControl.txt
  std::map<int, bool> nodeStructuralMap;  // boolean indicating structural BC
  std::map<int, bool> nodeMeshMotionMap;  // boolean indicating mesh motion BC
  std::map<int, bool> nodeThermalMap;     // boolean indicating heat transfer BC

  for (const auto &bcInfo : this->opts_.bcInfo) {
    int patchNum = bcInfo->patchNum;
    FaceBC *faceBC;
    NodeBC *nodeBC;
    if ((faceBC = dynamic_cast<FaceBC *>(bcInfo.get()))) {
      faceTypeMap[patchNum] = faceBC->rocFracFSIType;
    } else if ((nodeBC = dynamic_cast<NodeBC *>(bcInfo.get()))) {
      nodeTypeMap[patchNum] = nodeBC->rocfracControlType;
      nodeStructuralMap[patchNum] = nodeBC->structural;
      nodeMeshMotionMap[patchNum] = nodeBC->meshMotion;
      nodeThermalMap[patchNum] = nodeBC->thermal;
    }
  }

  // create meshBase object
  std::shared_ptr<meshBase> myMesh =
      meshBase::CreateShared(this->files_.inputMeshFile);

  // create PATRAN object from meshBase
  PATRAN::patran(myMesh, this->files_.inputMeshFile,
                 this->files_.outputMeshFile, faceTypeMap, nodeTypeMap,
                 nodeStructuralMap, nodeMeshMotionMap, nodeThermalMap,
                 this->opts_.nodePatchPreference);
  // pat->write();
}

}  // namespace DRV
}  // namespace NEM
