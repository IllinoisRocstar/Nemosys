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
#include "AuxiliaryFunctions.H"
#include "MeshPartitioning/meshStitcher.H"

#include "IO/cgnsAnalyzer.H"
#include "IO/rocstarCgns.H"
#include "Mesh/meshBase.H"

meshStitcher::meshStitcher(std::vector<std::string> _cgFileNames, bool surf)
    : cgFileNames(std::move(_cgFileNames)),
      stitchedMesh(nullptr),
      cgObj(nullptr) {
  if (!cgFileNames.empty()) {
    if (surf)
      initSurfCgObj();
    else
      initVolCgObj();
  }
}

void meshStitcher::initVolCgObj() {
  partitions.resize(cgFileNames.size(), nullptr);
  for (int iCg = 0; iCg < cgFileNames.size(); ++iCg) {
    partitions[iCg] = std::make_shared<cgnsAnalyzer>(cgFileNames[iCg]);
    partitions[iCg]->loadGrid(1);
    // defining partition flags
    std::vector<double> slnData(partitions[iCg]->getNElement(), iCg);
    // append partition number data if not already existing
    bool exists = false;
    std::vector<std::string> slnNamelist;
    partitions[iCg]->getSolutionDataNames(slnNamelist);
    for (const auto &in : slnNamelist)
      if (in == "partitionOld") {
        exists = true;
        break;
      }
    if (!exists)
      partitions[iCg]->appendSolutionData("partitionOld", slnData,
                                          solution_type_t::ELEMENTAL,
                                          partitions[iCg]->getNElement(), 1);
    // stitching new partitions to partition 0
    // close partition CGNS file to avoid clutter
    // MS: writing ghost mesh before closing
    // std::vector<std::string> secNames;
    // partitions[iCg]->getSectionNames(secNames);
    // auto its = std::find(secNames.begin(), secNames.end(), ":T4:virtual");
    // if (its != secNames.end())
    //  meshBase::Create(partitions[iCg]->getSectionMesh(*its),
    //                   "_virtual_" + nemAux::find_name(cgFileNames[iCg]) +
    //                   ".vtu")->write();
    // End of ghost mesh
    if (iCg) {
      partitions[0]->stitchMesh(partitions[iCg].get(), true);
      partitions[iCg]->closeCG();
    }
  }
  std::cout << "Meshes stitched successfully." << std::endl;
  std::cout << "Exporting stitched mesh to VTK format." << std::endl;
  std::string newname(cgFileNames[0]);
  std::size_t pos = newname.find_last_of('/');
  newname = newname.substr(pos + 1);
  newname = nemAux::trim_fname(newname, "stitched.vtu");
  stitchedMesh = meshBase::CreateShared(partitions[0]->getVTKMesh(), newname);
  std::cout << "Transferring physical quantities to vtk mesh." << std::endl;
  // figure out what is existing on the stitched grid
  int outNData, outNDim;
  std::vector<std::string> slnNameList;
  std::vector<std::string> appSlnNameList;
  partitions[0]->getSolutionDataNames(slnNameList);
  partitions[0]->getAppendedSolutionDataName(appSlnNameList);
  slnNameList.insert(slnNameList.end(), appSlnNameList.begin(),
                     appSlnNameList.end());

  // write all data into vtk file
  for (auto is = slnNameList.begin(); is < slnNameList.end(); is++) {
    std::vector<double> physData;
    partitions[0]->getSolutionDataStitched(*is, physData, outNData, outNDim);
    solution_type_t dt = partitions[0]->getSolutionDataObj(*is)->getDataType();
    if (dt == solution_type_t::NODAL) {
      std::cout << "Embedding nodal " << *is << std::endl;
      stitchedMesh->setPointDataArray((*is).c_str(), physData);
    } else {
      // gs field is 'weird' in irocstar files- we don't write it back
      // if (!(*is).compare("gs"))
      //  continue;
      std::cout << "Embedding cell-based " << *is << std::endl;
      stitchedMesh->setCellDataArray((*is).c_str(), physData);
    }
  }
  stitchedMesh->report();
  // TODO: Expensive IO, refactoring needed
  stitchedMesh->write();
}

void meshStitcher::initSurfCgObj() {
  cgObj = std::make_shared<rocstarCgns>(cgFileNames);
  // different read for burn files
  if (cgFileNames[0].find("burn") != std::string::npos)
    cgObj->setBurnBool(true);
  cgObj->loadCgSeries();
  cgObj->stitchGroup();
  // cgObj->closeCG();
  std::cout << "Surface mesh stitched successfully." << std::endl;
  std::cout << "Exporting stitched mesh to VTK format." << std::endl;
  std::string newname(cgFileNames[0]);
  std::size_t pos = newname.find_last_of('/');
  newname = newname.substr(pos + 1);
  newname = nemAux::trim_fname(newname, "stitched.vtu");
  stitchedMesh = meshBase::CreateShared(cgObj->getVTKMesh(), newname);
  std::cout << "Transferring physical quantities to vtk mesh." << std::endl;
  // figure out what exists on the stitched grid
  int outNData, outNDim;
  std::vector<std::string> slnNameList;
  std::vector<std::string> appSlnNameList;
  cgObj->getSolutionDataNames(slnNameList);
  cgObj->getAppendedSolutionDataName(appSlnNameList);
  slnNameList.insert(slnNameList.end(), appSlnNameList.begin(),
                     appSlnNameList.end());
  // write all data into vtk file
  for (auto is = slnNameList.begin(); is < slnNameList.end(); is++) {
    std::vector<double> physData;
    cgObj->getSolutionDataStitched(*is, physData, outNData, outNDim);
    solution_type_t dt = cgObj->getSolutionDataObj(*is)->getDataType();
    if (dt == solution_type_t::NODAL) {
      std::cout << "Embedding nodal " << *is << std::endl;
      stitchedMesh->setPointDataArray((*is).c_str(), physData);
    } else {
      // gs field is 'weird' in irocstar files- we don't write it back
      // if (!(*is).compare("gs"))
      if (*is == "mdot_old") { continue; }
      std::cout << "Embedding cell-based " << *is << std::endl;
      stitchedMesh->setCellDataArray((*is).c_str(), physData);
    }
  }
  stitchedMesh->report();
  // TODO: expensive IO, need to be refactored
  stitchedMesh->write();
}

std::shared_ptr<meshBase> meshStitcher::getStitchedMB() const {
  if (stitchedMesh)
    return stitchedMesh;
  else {
    std::cerr << "No stitched mesh to return!" << std::endl;
    exit(1);
  }
}

std::shared_ptr<cgnsAnalyzer> meshStitcher::getStitchedCGNS() const {
  if (!partitions.empty())
    return partitions[0];
  else if (cgObj)
    return cgObj;
  else {
    std::cerr << "No stitched mesh to return!" << std::endl;
    exit(1);
  }
}
