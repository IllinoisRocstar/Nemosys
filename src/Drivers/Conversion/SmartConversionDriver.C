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
#include "Drivers/Conversion/SmartConversionDriver.H"

#include "Mesh/geoMeshFactory.H"

namespace NEM {
namespace DRV {

SmartConversionDriver::SmartConversionDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

SmartConversionDriver::SmartConversionDriver()
    : SmartConversionDriver({{}, {}}, Opts{{}}) {}

SmartConversionDriver::Opts::Opts(jsoncons::optional<MeshData> meshData)
    : meshData(std::move(meshData)) {}

const SmartConversionDriver::Files &SmartConversionDriver::getFiles() const {
  return files_;
}

void SmartConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

void SmartConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

const SmartConversionDriver::Opts &SmartConversionDriver::getOpts() const {
  return opts_;
}

void SmartConversionDriver::execute() const {
  if (this->opts_.meshData.has_value()) {
    vtkSmartPointer<NEM::MSH::geoMeshBase> srcGM =
        NEM::MSH::Read(this->files_.inputMeshFile,
                       this->opts_.meshData.value().inputFormat.value_or(""));

    vtkSmartPointer<NEM::MSH::geoMeshBase> trgGM =
        NEM::MSH::New(this->files_.outputMeshFile,
                      this->opts_.meshData.value().outputFormat.value_or(""));

    trgGM->takeGeoMesh(srcGM);
    trgGM->write(this->files_.outputMeshFile);
  }

  else {
    vtkSmartPointer<NEM::MSH::geoMeshBase> srcGM =
        NEM::MSH::Read(this->files_.inputMeshFile);

    vtkSmartPointer<NEM::MSH::geoMeshBase> trgGM =
        NEM::MSH::New(this->files_.outputMeshFile);

    trgGM->takeGeoMesh(srcGM);
    trgGM->write(this->files_.outputMeshFile);
  }
}

}  // namespace DRV
}  // namespace NEM
