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
#include "Drivers/Conversion/GmshToVtkConversionDriver.H"

namespace NEM {
namespace DRV {

GmshToVtkConversionDriver::GmshToVtkConversionDriver(Files files)
    : files_(std::move(files)) {}

GmshToVtkConversionDriver::GmshToVtkConversionDriver()
    : GmshToVtkConversionDriver({{}, {}}) {}

const GmshToVtkConversionDriver::Files &GmshToVtkConversionDriver::getFiles()
    const {
  return files_;
}

void GmshToVtkConversionDriver::setFiles(Files files) {
  files_ = std::move(files);
}

void GmshToVtkConversionDriver::execute() const {
  if (this->files_.inputMeshFile.find(".msh") != std::string::npos) {
    std::cout << "Detected file in GMSH format" << std::endl;
    std::cout << "Converting to VTK ...." << std::endl;
  } else {
    std::cerr << "Source mesh file is not in GMSH format" << std::endl;
  }
  meshBase *mb = meshBase::exportGmshToVtk(this->files_.inputMeshFile);
  mb->report();
  mb->write(this->files_.outputMeshFile);
}

const GmshToVtkConversionDriver::Opts &GmshToVtkConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
