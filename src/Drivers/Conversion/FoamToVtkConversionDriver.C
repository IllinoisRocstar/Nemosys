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
#include "Drivers/Conversion/FoamToVtkConversionDriver.H"

#include "Mesh/foamMesh.H"
#include "Mesh/vtkMesh.H"

namespace NEM {
namespace DRV {

FoamToVtkConversionDriver::FoamToVtkConversionDriver(Files file)
    : file_(std::move(file)) {}

FoamToVtkConversionDriver::FoamToVtkConversionDriver()
    : FoamToVtkConversionDriver(Files{{}}) {}

const FoamToVtkConversionDriver::Files &FoamToVtkConversionDriver::getFiles()
const {
  return file_;
}

void FoamToVtkConversionDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

void FoamToVtkConversionDriver::execute() const {
  meshBase *fm = new FOAM::foamMesh();
  fm->read("NULL");
  // TODO: Fix report and write methods for the foamMesh class
  // std::cout << "Variable values is = " << srcmsh << std::endl;
  vtkMesh *vm = new vtkMesh(fm->getDataSet(), this->file_.outputFile);
  vm->report();
  vm->write();
  delete vm;
  delete fm;
}

const FoamToVtkConversionDriver::Opts &FoamToVtkConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
