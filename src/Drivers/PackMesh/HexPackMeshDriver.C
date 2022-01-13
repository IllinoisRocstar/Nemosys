/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#include "Drivers/PackMesh/HexPackMeshDriver.H"

#include <boost/filesystem.hpp>
#include "AuxiliaryFunctions.H"
#include "MeshManipulationFoam/MeshManipulationFoam.H"
#include "MeshGeneration/blockMeshGen.H"
#include "MeshGeneration/snappymeshGen.H"
#ifdef HAVE_GMSH
#include "Geometry/rocPack.H"
#endif

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
  auto objMsh = std::unique_ptr<MeshManipulationFoam>(
      new MeshManipulationFoam(&mmfParams));
  const char *nameFile = "a";  // Dummy name for input

  // A rocpack method that creates stl and then moves it to triSurface using
  // boost.
  if (this->files_.isInputRocpackFile()) {
#ifdef HAVE_GMSH
    std::string hexOutSTL = this->files_.getInputFile() + ".stl";
    auto objrocPck = std::unique_ptr<NEM::GEO::rocPack>(
        new NEM::GEO::rocPack(this->files_.getInputFile(), hexOutSTL));

    objrocPck->removeBoundaryVolumes();
    objrocPck->rocPack2Surf();

    // const std::string dir_path11 = hexOutSTL;
    boost::filesystem::path dir11(hexOutSTL);

    const std::string dir_path2 =
        "./constant/triSurface/" + this->files_.getInputFile() + ".stl";
    boost::filesystem::path dir2(dir_path2);

    boost::filesystem::copy_file(
        dir11, dir2, boost::filesystem::copy_option::overwrite_if_exists);
#else
    std::cerr << "Cannot process rocPack input file without gmsh.\n";
    std::exit(1);
#endif
  } else {
    const std::string dir_path11 = this->files_.getInputFile();
    boost::filesystem::path dir11(dir_path11);

    const std::string dir_path2 =
        "./constant/triSurface/" + this->files_.getInputFile();
    boost::filesystem::path dir2(dir_path2);

    boost::filesystem::copy_file(
        dir11, dir2, boost::filesystem::copy_option::overwrite_if_exists);
  }

  // blockMesh utility takes user input for surrounding box region and
  // generates mesh block in constant/polyMesh folder. It will overwrite the
  // previous mesh created by CfMesh. This mesh will be used as background mesh
  // by snappyHexMesh later.
  auto bmParamsCopy = this->opts_.bmParams;  // should outlive objBM
  bmParamsCopy.isPackMesh = true;
  auto objBM = std::unique_ptr<blockMeshGen>(new blockMeshGen(&bmParamsCopy));
  objBM->createMeshFromSTL(nameFile);

  // Updating location for next process
  auto adjust = this->opts_.locAdjust >= 0. ? this->opts_.locAdjust : 0.;
  auto smParamsCopy = this->opts_.smParams;  // should outlive objSHM
  smParamsCopy.isPackMesh = true;
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
  auto objSHM =
      std::unique_ptr<snappymeshGen>(new snappymeshGen(&smParamsCopy));
  objSHM->createMeshFromSTL(nameFile);

  // splitMeshRegions reads mesh from constant/polyMesh directory and splits
  // mesh into separate regions. Each region will represent one solid. It will
  // write bunch of domain.* directories inside constant/ and system/ folders.
  // This function does a random walking around mesh to identify different
  // regions and in that process, while naming all domains as domain.*
  // in constant folder, it skips one number for the disconnected region it
  // encounters first. This number is taken out to provide as input to merge
  // mesh.
  std::pair<std::vector<int>, std::vector<std::string>> dirStat =
      objMsh->splitMshRegions();

  std::string surroundingRegion = dirStat.second[0];
  mmfParams.surfSplitParams.pckRegionNames = dirStat.second;

  /*std::pair<std::vector<int>, std::string> dirStat =
  objMsh->splitMshRegions(); int skippedDir = dirStat.first[0]; int totalRegs =
  dirStat.first[1] - 1; std::string surroundingRegion = dirStat.second;
  */

  // *************** Merge meshes using geoMeshBase and write
  // *******************// packRegNames is a vector containing all names of pack
  // regions
  mmfParams.surfSplitParams.pckRegionNames.erase(
      mmfParams.surfSplitParams.pckRegionNames.begin());
  std::vector<std::string> packRegNames =
      mmfParams.surfSplitParams.pckRegionNames;

  // FoamGeoMesh
  auto fgm_ = NEM::MSH::Read(surroundingRegion + ".foam");
  auto mesh = NEM::MSH::New(this->files_.outCombinedFile);
  mesh->takeGeoMesh(fgm_);

  for (auto &packRegName : packRegNames) {
    auto fgm_loop = NEM::MSH::Read(packRegName + ".foam");
    mesh->mergeGeoMesh(fgm_loop);
  }

  mesh->write(this->files_.outCombinedFile);

  mesh->Delete();
  fgm_->Delete();

  // *************** Merge meshes using meshBase and write *******************//

  //  bool readDB = false;
  //  // Create surrounding region database
  //  meshBase *fm = new FOAM::foamMesh(readDB);
  //  fm->read(surroundingRegion);
  //  std::vector<double> physId = std::vector<double>(fm->getNumberOfCells(),
  //  0); auto *vm = new vtkMesh(fm->getDataSet(),
  //  this->files_.outCombinedFile); vm->setCellDataArray("PhysGrpId", physId);
  //
  //  // Loop through all pack particles and merge their databases into main
  //  // database
  //  for (int i = 0; i < (int) packRegNames.size(); i++) {
  //    fm->read(packRegNames[i]);
  //    std::vector<double> physIdLoop =
  //        std::vector<double>(fm->getNumberOfCells(), i + 1);
  //    auto *vm2 =
  //        new vtkMesh(fm->getDataSet(), "pack" + std::to_string(i) + ".vtu");
  //    vm2->setCellDataArray("PhysGrpId", physIdLoop);
  //    // vm2->write();
  //    vm->merge(vm2->getDataSet());
  //  }
  //
  //  // Write mesh and clean up objects
  //  vm->report();
  //  vm->write();
  //  delete vm;
  //  delete fm;

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

  // End of workflow
}

}  // namespace DRV
}  // namespace NEM
