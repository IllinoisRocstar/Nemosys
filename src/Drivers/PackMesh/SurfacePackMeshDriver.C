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
#include "Drivers/PackMesh/SurfacePackMeshDriver.H"

#include "Geometry/rocPack.H"

namespace NEM {
namespace DRV {

SurfacePackMeshDriver::Files::Files(std::string rocpackFile,
                                    std::string outputMeshFile)
    : rocpackFile(std::move(rocpackFile)),
      outputMeshFile(std::move(outputMeshFile)) {}

SurfacePackMeshDriver::Opts::Opts(
    jsoncons::optional<Periodic3DOpts> periodic3DOpts)
    : periodic3DOpts(std::move(periodic3DOpts)) {}

SurfacePackMeshDriver::SurfacePackMeshDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

SurfacePackMeshDriver::SurfacePackMeshDriver()
    : SurfacePackMeshDriver({{}, {}}, Opts{{}}) {}

const SurfacePackMeshDriver::Files &SurfacePackMeshDriver::getFiles() const {
  return files_;
}

void SurfacePackMeshDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const SurfacePackMeshDriver::Opts &SurfacePackMeshDriver::getOpts() const {
  return opts_;
}

void SurfacePackMeshDriver::setOpts(Opts opts) {
  if (opts.periodic3DOpts.has_value()) {
    const auto &periodic3DOpts = opts.periodic3DOpts.value();
    if (!(opts.meshAlgorithm == 1 || opts.meshAlgorithm == 4 ||
          opts.meshAlgorithm == 7 || opts.meshAlgorithm == 9 ||
          opts.meshAlgorithm == 10)) {
      std::cerr << "Valid choices for 3D meshing algorithm are "
                << " 1 (Delunay), 4 (Frontal), 7 (MMG3D), 9 (R-Tree), 10 (HXT)"
                << std::endl;
      throw;
    }
    if (!(periodic3DOpts.elemOrder == 1 || periodic3DOpts.elemOrder == 2)) {
      std::cerr << "Only element orders 1 and 2 are supported!" << std::endl;
      throw;
    }
  } else {
    if (opts.meshAlgorithm == 1 || opts.meshAlgorithm == 2 ||
        opts.meshAlgorithm == 5 || opts.meshAlgorithm == 6 ||
        opts.meshAlgorithm == 7 || opts.meshAlgorithm == 8 ||
        opts.meshAlgorithm == 9) {
    } else {
      std::cerr
          << "Valid choices for 2D meshing algorithm are "
          << " 1 (MeshAdapt), 2 (Automatic),  5 (Delunay), "
          << " 6 (Frontal Delunay), 7 (BAMG), 8 (Frontal Delunay for Quads)"
          << " 9 (Packing of Parallelograms)." << std::endl;
      throw;
    }
  }
  this->opts_ = std::move(opts);
}

// TODO
// 1. Some implementation to improve mesh if mesh is not of desired
//    quality.
void SurfacePackMeshDriver::execute() const {
  auto *objrocPck = new NEM::GEO::rocPack(this->files_.rocpackFile,
                                          this->files_.outputMeshFile);

  if (this->opts_.upperThreshold.has_value() ||
      this->opts_.lowerThreshold.has_value())
    objrocPck->applyFilter(this->opts_.upperThreshold.value_or(0.),
                           this->opts_.lowerThreshold.value_or(0.));

  objrocPck->setMeshingAlgorithm(this->opts_.meshAlgorithm);

  if (this->opts_.preserveSize) objrocPck->setSizePreservation();

  if (this->opts_.refineLevel.has_value())
    objrocPck->assignRefinement(this->opts_.refineLevel.value());

  if (this->opts_.removeBoundaryPacks) objrocPck->removeBoundaryVolumes();

  if (this->opts_.scaleValue.has_value())
    objrocPck->shrinkVolumes(this->opts_.scaleValue.value());

  if (this->opts_.meshSize.has_value())
    objrocPck->setMeshSize(this->opts_.meshSize.value());

  if (this->opts_.enableDefaultOut) objrocPck->enableDefOuts();

  // Create periodic mesh
  if (this->opts_.periodic3DOpts.has_value()) {
    const auto &periodic3D = this->opts_.periodic3DOpts.value();
    if (periodic3D.customDomain && periodic3D.setPeriodicGeo) {
      std::cerr
          << "WARNING!! -> Cannot make geometry periodic using custom bounds."
          << " using remove geometries on boundary option instead for meshing!"
          << std::endl;
    }
    objrocPck->translateAll(periodic3D.transferMesh[0],
                            periodic3D.transferMesh[1],
                            periodic3D.transferMesh[2]);

    objrocPck->setElementOrder(periodic3D.elemOrder);

    if (periodic3D.createCohesive) {
      objrocPck->sanityCheckOn();
      objrocPck->enableCohesiveElements();
    }

    switch (periodic3D.physGrpOptions) {
      case PhysGrpOpts::MULTI: objrocPck->enablePhysicalGrps(); break;
      case PhysGrpOpts::TWO: objrocPck->enableTwoPhysGrps(); break;
      case PhysGrpOpts::PER_SHAPE:
        objrocPck->enablePhysicalGroupsPerShape();
        break;
      case PhysGrpOpts::NONE: break;
    }

    if (periodic3D.enablePatches) objrocPck->enableSurfacePatches();

    if (periodic3D.setPeriodicGeo) objrocPck->setPeriodicGeometry();

    if (periodic3D.customDomain) {
      const auto &customDomain = periodic3D.customDomain.value();
      objrocPck->setCustomDomain(
          {customDomain.initial[0], customDomain.initial[1],
           customDomain.initial[2], customDomain.length[0],
           customDomain.length[1], customDomain.length[2]});
      objrocPck->removeBoundaryVolumes();
      objrocPck->setPeriodicMesh();
      objrocPck->rocPack2Periodic3D();
    } else {
      objrocPck->setPeriodicGeometry();
      objrocPck->setPeriodicMesh();
      objrocPck->rocPack2Periodic3D();
    }

  } else {  // Generate STL only
    objrocPck->rocPack2Surf();
  }

  if (objrocPck) delete objrocPck;
}

}  // namespace DRV
}  // namespace NEM
