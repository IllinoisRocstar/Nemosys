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
#include "Drivers/NucMeshDriver.H"

#include "Mesh/geoMeshFactory.H"
#include "Mesh/smeshGeoMesh.H"
#include "Services/NucMeshSrv.H"

namespace NEM {
namespace DRV {

NucMeshDriver::NucMeshDriver(Files file, Opts opts)
    : file_(std::move(file)), opts_(std::move(opts)) {}

NucMeshDriver::NucMeshDriver() : NucMeshDriver(Files{std::string{}}, Opts{}) {}

const NucMeshDriver::Files &NucMeshDriver::getFiles() const { return file_; }

void NucMeshDriver::setFiles(Files files) { file_ = std::move(files); }

const NucMeshDriver::Opts &NucMeshDriver::getOpts() const { return opts_; }

void NucMeshDriver::setOpts(Opts opts) { opts_ = std::move(opts); }

jsoncons::string_view NucMeshDriver::getProgramType() const {
  return programType;
}

vtkSmartPointer<NEM::MSH::geoMeshBase> NucMeshDriver::draw() const {
  vtkNew<NEM::SRV::NucMeshSrv> nucMeshRunner{};
  nucMeshRunner->SetConfiguration(opts_);
  nucMeshRunner->Update();
  return NEM::MSH::smeshGeoMesh::SafeDownCast(nucMeshRunner->GetOutput());
}

void NucMeshDriver::execute() const {
  auto outNative = this->draw();
  auto outType = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::New(file_.outputFile));
  outType->takeGeoMesh(outNative);
  outType->write(file_.outputFile);
}

}  // namespace DRV
}  // namespace NEM