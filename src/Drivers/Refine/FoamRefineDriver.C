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
#include "Drivers/Refine/FoamRefineDriver.H"

#include <pointFieldsFwd.H>
#include <interpolatePointToCell.H>
#include "getDicts.H"
#include "Refinement/AMRFoam.H"
#include "AuxiliaryFunctions.H"
#include "Mesh/foamGeoMesh.H"
#include "Mesh/geoMeshFactory.H"

#ifdef MLAMR
#  include <Mesh/meshBase.H>
#  include <fdeep/fdeep.hpp>
#endif

// helpers for execute
namespace {
std::unique_ptr<Foam::AMRFoam> initializeRefine(
    const NEM::DRV::RefineDriver::Files &files,
    const NEM::DRV::FoamRefineOptsBase &opts, const Foam::Time &runTime) {
  // Reading mesh from vtk file.
  auto vgm = NEM::MSH::Read(files.inputMeshFile);
  auto mesh = new NEM::MSH::foamGeoMesh();
  mesh->takeGeoMesh(vgm);
  mesh->write("");
  mesh->Delete();
  vgm->Delete();

  Foam::polyMesh mesh1(Foam::IOobject(Foam::polyMesh::defaultRegion,
                                      runTime.timeName(), runTime,
                                      Foam::IOobject::MUST_READ));

  // Initialized AMR class
  auto amr = std::unique_ptr<Foam::AMRFoam>{new Foam::AMRFoam(mesh1)};

  if (opts.writeFieldData)
    if (!opts.writeMesh) amr->enableUpdatedField();

  if (opts.writeRefHistory)
    if (!opts.writeMesh) amr->enableRefHistoryData();

  if (opts.writeMesh) amr->enableMeshWriting();

  return amr;
}

void finishRefine(Foam::AMRFoam *amr, const std::string &outputMeshFile,
                  const Foam::Time &runTime) {
  amr->writeMesh();

  // Export final mesh to vtu format
  auto fm = new Foam::fvMesh(Foam::IOobject("", runTime.timeName(), runTime,
                                            Foam::IOobject::MUST_READ));
  auto fgm_ = new NEM::MSH::foamGeoMesh(fm);
  auto mesh = NEM::MSH::New(outputMeshFile);
  mesh->takeGeoMesh(fgm_);
  mesh->write(outputMeshFile);
  mesh->Delete();
  fgm_->Delete();
}
}  // namespace

namespace NEM {
namespace DRV {

FoamRefineDriver::Opts::Opts(Criteria refCriteria, std::string inputFieldFile)
    : FoamRefineOptsBase(),
      inputFieldFile(std::move(inputFieldFile)),
      refCriteria(refCriteria) {}

FoamRefineDriver::FoamRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

FoamRefineDriver::FoamRefineDriver() : FoamRefineDriver({{}, {}}, {{}, {}}) {}

const FoamRefineDriver::Opts &FoamRefineDriver::getOpts() const {
  return opts_;
}

void FoamRefineDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

void FoamRefineDriver::execute() const {
  // Initializes AMR workflow
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  initFoam->writeBasicDicts("system/", this->opts_.startTime,
                            this->opts_.timeStep, this->opts_.endTime);
  initFoam->createDynamicMeshDict(true);
  auto controlDict_ = initFoam->createControlDict(true);
  Foam::Time runTime(controlDict_.get(), ".", ".");

  auto amr = initializeRefine(this->files_, this->opts_, runTime);

  int nSteps = static_cast<int>(std::round(opts_.endTime / opts_.timeStep));

  // Creates scalar field from incoming text/CSV files
  Foam::volScalarField meshFieldXY =
      amr->readIncomingCellField(this->opts_.inputFieldFile, "meshFieldXY");

  if (opts_.refCriteria == Opts::Criteria::VALUE) {
    for (int i = 0; i < nSteps; i++) {
      amr->updateAMR(this->opts_.refineInterval, this->opts_.maxRefinement,
                     meshFieldXY, this->opts_.lowerRefineLevel,
                     this->opts_.upperRefineLevel, this->opts_.unrefineAbove,
                     this->opts_.unrefineBelow, this->opts_.nBufferLayers,
                     this->opts_.maxCells);

      runTime++;
    }
    meshFieldXY.write();
  } else {  // opts_.refCriteria == Opts::Criteria::GRADIENT
    Foam::volScalarField magGrad = amr->getGradient(meshFieldXY);
    for (int i = 0; i < nSteps; i++) {
      amr->updateAMR(this->opts_.refineInterval, this->opts_.maxRefinement,
                     magGrad, this->opts_.lowerRefineLevel,
                     this->opts_.upperRefineLevel, this->opts_.unrefineAbove,
                     this->opts_.unrefineBelow, this->opts_.nBufferLayers,
                     this->opts_.maxCells);

      runTime++;
    }
    magGrad.write();
  }

  finishRefine(amr.get(), this->files_.outputMeshFile, runTime);

  std::cout << "End!" << std::endl;
}

#ifdef MLAMR
FoamMLRefineDriver::FoamMLRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

FoamMLRefineDriver::FoamMLRefineDriver() : FoamMLRefineDriver({{}, {}}, {}) {}

const FoamMLRefineDriver::Opts &FoamMLRefineDriver::getOpts() const {
  return opts_;
}

void FoamMLRefineDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

void FoamMLRefineDriver::execute() const {
  // Initializes AMR workflow
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  initFoam->writeBasicDicts("system/", this->opts_.startTime,
                            this->opts_.timeStep, this->opts_.endTime);
  initFoam->createDynamicMeshDict(true);
  auto controlDict_ = initFoam->createControlDict(true);
  Foam::Time runTime(controlDict_.get(), ".", ".");

  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  auto amr = initializeRefine(this->files_, this->opts_, runTime);

  int nSteps = opts_.endTime / opts_.timeStep;

  std::vector<double> nonDimUGrad;
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  mesh->getCellDataArray("nonDimUGrad", nonDimUGrad);
  mesh->getCellDataArray("X", X);
  mesh->getCellDataArray("Y", Y);
  mesh->getCellDataArray("Z", Z);

  // Loading ML Model
  const auto model = fdeep::load_model(this->opts_.mlModel);

  // Making predictions
  std::vector<int> refinementVec;
  for (int i = 0; i < nonDimUGrad.size(); i++) { refinementVec.push_back(0); }
  for (int i = 0; i < nonDimUGrad.size(); i++) {
    const auto result = model.predict(
        {fdeep::tensor(fdeep::tensor_shape(static_cast<double>(4)),
                       {X[i], Y[i], Z[i], nonDimUGrad[i]})});

    if (result[0].get(0, 0, 0, 0, 0) >= 0.5) refinementVec[i] = 1;
  }
  Foam::volScalarField meshField = amr->assignToVolScalarField(refinementVec);

  for (int i = 0; i < nSteps; i++) {
    amr->updateAMRML(this->opts_.refineInterval, this->opts_.maxRefinement,
                     this->opts_.nBufferLayers, this->opts_.maxCells,
                     meshField);

    runTime++;
  }

  finishRefine(amr.get(), this->files_.outputMeshFile, runTime);
  std::cout << "End!" << std::endl;
}
#endif

}  // namespace DRV
}  // namespace NEM
