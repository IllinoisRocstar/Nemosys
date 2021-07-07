#include "Drivers/PackMesh/HexPackMeshDriver.H"

#include <boost/filesystem.hpp>
#include "AuxiliaryFunctions.H"
#include "MeshManipulationFoam/MeshManipulationFoam.H"
#include "MeshGeneration/blockMeshGen.H"
#include "Geometry/rocPack.H"
#include "MeshGeneration/snappymeshGen.H"
#include "Mesh/vtkMesh.H"

namespace NEM {
namespace DRV {

HexPackMeshDriver::Files::Files(std::string inputFile, bool isRocpack,
                                std::string outPackMeshFile,
                                std::string outSurroundingFile)
    : outPackMeshFile(std::move(outPackMeshFile)),
      outSurroundingFile(std::move(outSurroundingFile)),
      rocpackOrGeoFile(std::move(inputFile)),
      useRocpack(isRocpack) {}

void HexPackMeshDriver::Files::setRocpackFile(std::string rocpackFile) {
  rocpackOrGeoFile = std::move(rocpackFile);
  useRocpack = true;
}

void HexPackMeshDriver::Files::setGeoFile(std::string geoFile) {
  rocpackOrGeoFile = std::move(geoFile);
  useRocpack = false;
}

const std::string &HexPackMeshDriver::Files::getInputFile() const {
  return rocpackOrGeoFile;
}

bool HexPackMeshDriver::Files::isInputRocpackFile() const { return useRocpack; }

HexPackMeshDriver::HexPackMeshDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

HexPackMeshDriver::HexPackMeshDriver()
    : HexPackMeshDriver({{}, {}, {}, {}}, {}) {}

const HexPackMeshDriver::Files &HexPackMeshDriver::getFiles() const {
  return files_;
}

void HexPackMeshDriver::setFiles(Files files) {
  const auto stlFile = files.isInputRocpackFile()
                           ? files.getInputFile() + ".stl"
                           : files.getInputFile();
  this->opts_.smParams.geomFileName = stlFile;
  if (auto box = dynamic_cast<bmBox *>(this->opts_.bmParams.shape.get())) {
    if (box->autoGenerate.has_value()) {
      box->autoGenerate.value().packFileName = stlFile;
    }
  }
  this->files_ = std::move(files);
}

const HexPackMeshDriver::Opts &HexPackMeshDriver::getOpts() const {
  return opts_;
}

void HexPackMeshDriver::setOpts(Opts opts) {
  const auto stlFile = getFiles().isInputRocpackFile()
                           ? getFiles().getInputFile() + ".stl"
                           : getFiles().getInputFile();
  if (auto box = dynamic_cast<bmBox *>(opts.bmParams.shape.get())) {
    if (box->autoGenerate.has_value()) {
      if (!box->autoGenerate.value().packFileName.empty()) {
        std::cerr
            << "opts.bmParams.shape.autoGenerate.packFileName is ignored. "
               "Using getFiles().getInputFile() for geometry file"
            << std::endl;
      }
      box->autoGenerate.value().packFileName = stlFile;
    }
  } else {
    std::cerr << "HexPackMeshDriver only accepts box shape for blockmesh.\n";
    exit(1);
  }

  if (!opts.smParams.geomFileName.empty()) {
    std::cerr << "opts.smParams.geomFileName is ignored. Using "
                 "getFiles().getInputFile() for geometry file"
              << std::endl;
  }
  opts.smParams.geomFileName = stlFile;

  // Default values specific to this driver
  if (opts.mmfMergeParams.masterCasePath.empty()) {
    opts.mmfMergeParams.masterCasePath = ".";
  }
  if (opts.mmfMergeParams.addCasePath.empty()) {
    opts.mmfMergeParams.addCasePath = ".";
  }
  if (opts.mmfCreatePatchParams.surroundingName.empty()) {
    opts.mmfCreatePatchParams.surroundingName = "Soil";
  }
  if (opts.mmfCreatePatchParams.packsName.empty()) {
    opts.mmfCreatePatchParams.packsName = "Rocks";
  }
  if (opts.mmfCreatePatchParams.srrndngPatchType.empty()) {
    opts.mmfCreatePatchParams.srrndngPatchType = "wall";
  }
  if (opts.mmfCreatePatchParams.packsPatchType.empty()) {
    opts.mmfCreatePatchParams.packsPatchType = "wall";
  }
  this->opts_ = std::move(opts);
}

void HexPackMeshDriver::execute() const {
  // Makes sure that constant and triSurface directories are present.
  // If directories are already present, it will not do anything.
  const char dir_path[] = "./constant";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating triSurface directory"
              << "for the snappyHexMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  const char dir_path1[] = "./constant/triSurface";
  boost::filesystem::path dir1(dir_path1);
  try {
    boost::filesystem::create_directory(dir1);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating triSurface directory"
              << "for the snappyHexMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  // This object has access to all MeshManipulation utilities.
  MeshManipulationFoamParams mmfParams;  // should outlive objMsh
  mmfParams.mergeMeshesParams = this->opts_.mmfMergeParams;
  mmfParams.createPatchParams = this->opts_.mmfCreatePatchParams;
  MeshManipulationFoam *objMsh = new MeshManipulationFoam(&mmfParams);
  const char *nameFile = "a";  // Dummy name for input

  // A rocpack method that creates stl and then moves it to triSurface using
  // boost.
  if (this->files_.isInputRocpackFile()) {
    std::string hexOutSTL = this->files_.getInputFile() + ".stl";
    auto *objrocPck =
        new NEM::GEO::rocPack(this->files_.getInputFile(), hexOutSTL);

    objrocPck->removeBoundaryVolumes();
    objrocPck->rocPack2Surf();

    if (objrocPck) delete objrocPck;

    const std::string dir_path1 = hexOutSTL;
    boost::filesystem::path dir1(dir_path1);

    const std::string dir_path2 =
        "./constant/triSurface/" + this->files_.getInputFile() + ".stl";
    boost::filesystem::path dir2(dir_path2);

    boost::filesystem::copy_file(
        dir1, dir2, boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    const std::string dir_path1 = this->files_.getInputFile();
    boost::filesystem::path dir1(dir_path1);

    const std::string dir_path2 =
        "./constant/triSurface/" + this->files_.getInputFile();
    boost::filesystem::path dir2(dir_path2);

    boost::filesystem::copy_file(
        dir1, dir2, boost::filesystem::copy_option::overwrite_if_exists);
  }

  // blockMesh utility takes user input for surrounding box region and
  // generates mesh block in constant/polyMesh folder. It will overwrite the
  // previous mesh created by CfMesh. This mesh will be used as background mesh
  // by snappyHexMesh later.
  auto bmParamsCopy = this->opts_.bmParams;  // should outlive objBM
  blockMeshGen *objBM = new blockMeshGen(&bmParamsCopy);
  objBM->createMeshFromSTL(nameFile);

  // Updating location for next process
  auto adjust = this->opts_.locAdjust >= 0. ? this->opts_.locAdjust : 0.;
  auto smParamsCopy = this->opts_.smParams;  // should outlive objSHM
  smParamsCopy.castMeshControls.locMesh[0] =
      std::dynamic_pointer_cast<bmBox>(bmParamsCopy.shape)->coordsBox.first[0] +
      0.001 + adjust;
  smParamsCopy.castMeshControls.locMesh[1] =
      std::dynamic_pointer_cast<bmBox>(bmParamsCopy.shape)->coordsBox.first[1] +
      0.001 + adjust;
  smParamsCopy.castMeshControls.locMesh[2] =
      std::dynamic_pointer_cast<bmBox>(bmParamsCopy.shape)->coordsBox.first[2] +
      0.001 + adjust;

  // snappyHexMesh reads background mesh created using blockMesh and takes
  // surface file from foamToSurface utility to snap pack surface onto back-
  // ground mesh and creates different cellZones (i.e different solids). These
  // interfaces between pack and surrounding regions are completely conformal
  // due to snappyHexMesh's unique snapping abilities.
  snappymeshGen *objSHM = new snappymeshGen(&smParamsCopy);
  objSHM->createMeshFromSTL(nameFile);

  // splitMeshRegions reads mesh from constant/polyMesh directory and splits
  // mesh into separate regions. Each region will represent one solid. It will
  // write bunch of domain.* directories inside constant/ and system/ folders.
  // This function does a random walking around mesh to identify different
  // regions and in that process, while naming all domains as domain.*
  // in constant folder, it skips one number for the disconnected region it
  // encounters first. This number is taken out to provide as input to merge
  // mesh.
  std::pair<std::vector<int>, std::string> dirStat = objMsh->splitMshRegions();

  int skippedDir = dirStat.first[0];
  int totalRegs = dirStat.first[1] - 1;
  std::string surroundingRegion = dirStat.second;

  // *************** Merge meshes using meshBase and write *******************//
  // packRegNames is a vector containing all names of pack regions
  std::vector<std::string> packRegNames =
      mmfParams.surfSplitParams.pckRegionNames;

  bool readDB = false;
  // Create surrounding region database
  meshBase *fm = new FOAM::foamMesh(readDB);
  fm->read(surroundingRegion);
  std::vector<double> physId = std::vector<double>(fm->getNumberOfCells(), 0);
  vtkMesh *vm = new vtkMesh(fm->getDataSet(), this->files_.outCombinedFile);
  vm->setCellDataArray("PhysGrpId", physId);

  // Loop through all pack particles and merge their databases into main
  // database
  for (int i = 0; i < packRegNames.size(); i++) {
    fm->read(packRegNames[i]);
    std::vector<double> physIdLoop =
        std::vector<double>(fm->getNumberOfCells(), i + 1);
    vtkMesh *vm2 = new vtkMesh(fm->getDataSet(), this->files_.outCombinedFile);
    vm2->setCellDataArray("PhysGrpId", physIdLoop);
    vm2->write();
    vm->merge(vm2->getDataSet());
  }

  // Write mesh and clean up objects
  vm->write();
  if (vm) delete vm;
  if (fm) delete fm;
  // *************** Merge meshes using meshBase and write *******************//

  // *************** Merge meshes using OF and write ************************ //
  // Use this method for writing pack and surrounding meshes separately. This
  // uses mergeMeshes method from "src/MeshManipulation/MeshManipulationFoam.C".
  // It also writes merged mesh at the end. Use this as a replacement for the
  // method directly above (Merge meshes using meshBase).

  // mergeMeshes will read master domain (defined by user) from constant folder
  // and start merging other domain to it untill all slave domains are attached
  // to master domain. It loops through all domains in sequential manner and
  // skips the missing domain.* (also skipped by splitMeshRegions) to avoid
  // runtime error.
  // std::cout << "Total # of domains are = " << totalRegs << std::endl;
  // if (totalRegs == 1) {
  //   // Nothing
  // } else {
  //   objMsh->mergeMeshes(skippedDir, totalRegs);
  // }

  // // Reads current mesh and write it to separate VTK/VTU files
  // // Reads and converts pack mesh
  // bool readDB = false;
  // std::string regNme;
  // if (skippedDir == 1)
  //   regNme = "domain2";
  // else
  //   regNme = "domain1";

  // if (totalRegs == 1) regNme = _snappyparams->singleSolidPatch;

  // meshBase *fm = new FOAM::foamMesh(readDB);
  // fm->read(regNme);
  // vtkMesh *vm2 = new vtkMesh(fm->getDataSet(), ofname_pack);
  // std::vector<double> physIdPack
  //   = std::vector<double>(fm->getNumberOfCells(),1);
  // vm2->setCellDataArray("PhysGrpId",physIdPack);
  // vm2->write();

  // // Reads and converts surronding mesh
  // regNme = surroundingRegion;

  // if (totalRegs == 1) regNme = surroundingRegion;
  // meshBase *fm2 = new FOAM::foamMesh(readDB);
  // fm2->read(regNme);
  // std::vector<double> physIdSurrounding
  //   = std::vector<double>(fm2->getNumberOfCells(),0);
  // vtkMesh *vm3 = new vtkMesh(fm2->getDataSet(), ofname_surrndng);
  // vm3->setCellDataArray("PhysGrpId",physIdSurrounding);
  // vm3->write();

  // // Merge meshes and write
  // vtkMesh *vm = new vtkMesh(vm2->getDataSet(),ofname_merged);
  // vm->merge(vm3->getDataSet());
  // vm->report();
  // vm->write();
  // if (vm) delete vm;
  // if (fm) delete fm;
  // if (vm2) delete vm2;
  // if (fm2) delete fm2;
  // if (vm3) delete vm3;

  // // Outputs useful mesh quality parameters for users
  // std::string SurroundingName = "surroundingMeshQuality";
  // std::string Surroundingmesh = "geom_surrounding_mesh.vtu";
  // std::string PackName = "packMeshQuality";
  // std::string Packmesh = "geom_pack_mesh.vtu";
  // MeshQualityDriver *objSurrQ =
  //     new MeshQualityDriver(Surroundingmesh, SurroundingName);
  // MeshQualityDriver *objPackQ = new MeshQualityDriver(Packmesh, PackName);
  // if (objSurrQ) delete objSurrQ;
  // if (objPackQ) delete objPackQ;
  // *************** Merge meshes using OF and write ************************ //

  // Cleaning up
  if (objMsh) delete objMsh;
  if (objSHM) delete objSHM;
  if (objBM) delete objBM;
  // End of workflow
}

}  // namespace DRV
}  // namespace NEM
