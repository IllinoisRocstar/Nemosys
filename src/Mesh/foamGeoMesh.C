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
#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/foamGeoMesh.H"
#include "getDicts.H"

#include <iostream>
#include <set>
#include <string>
#include <utility>

#ifdef HAVE_GMSH
#  include <gmsh.h>
#endif
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkStringArray.h>

#include <IOobjectList.H>
#include <cellZone.H>
#include <cellZoneSet.H>
#include <fileName.H>
#include <foamVtkVtuAdaptor.H>
#include <globalMeshData.H>
#include <timeSelector.H>
#include <topoSetSource.H>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace MSH {

vtkStandardNewMacro(foamGeoMesh)

foamGeoMesh *foamGeoMesh::Read(const std::string &fileName,
                               const std::string &geoArrayName) {
  auto foamGM = new foamGeoMesh();
  Foam::word regionName;

  if (fileName.empty()) {
    regionName = Foam::fvMesh::defaultRegion;
  } else {
    regionName = fileName;
  }

  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  foamGM->controlDict_ = initFoam->createControlDict(writeDicts);

  foamGM->runTime_ = std::unique_ptr<Foam::Time>(
      new Foam::Time(foamGM->controlDict_.get(), ".", "."));

  foamGM->fmesh_.reset(new Foam::fvMesh(
      Foam::IOobject(regionName, foamGM->runTime_->timeName(),
                     *foamGM->runTime_, Foam::IOobject::MUST_READ)));

  auto gmMesh = foam2GM(foamGM->fmesh_.get());

  foamGM->setGeoMesh(gmMesh);

  return foamGM;
}

foamGeoMesh::foamGeoMesh() : foamGeoMesh(nullptr) { InitializeFoam(); }

foamGeoMesh::foamGeoMesh(Foam::fvMesh *foamMesh,
                         const std::string &phyGrpArrayName)
    : geoMeshBase(foam2GM(foamMesh, phyGrpArrayName)), fmesh_(foamMesh) {
  std::cout << "foamGeoMesh Constructed" << std::endl;
}

foamGeoMesh::~foamGeoMesh() {
  std::cout << "foamGeoMesh destructed" << std::endl;
}

geoMeshBase::GeoMesh foamGeoMesh::foam2GM(Foam::fvMesh *foamMesh,
                                          const std::string &phyGrpArrayName) {
  // create vtk database
  vtkSmartPointer<vtkUnstructuredGrid> vtkdataSet =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Check for cellsets and cellzones
  if (!foamMesh) {
    return {vtkdataSet, "", "", {}};
  } else {
    const Foam::cellZoneMesh &cellZones =
        foamMesh->cellZones();  // get cellZoneMesh
    Foam::label zoneI = cellZones.size();

    // creating equivalent vtk topology from fvMesh
    // by default polyhedral cells will be decomposed to
    // tets and pyramids. Additional points will be added
    // to underlying fvMesh.
    std::cout << "Performing topological decomposition.\n";
    auto objVfoam = Foam::vtk::vtuAdaptor();
    vtkdataSet = objVfoam.internal(*foamMesh);

    // Check for cellZones and extract physical groups
    Foam::List<Foam::scalar> physGrps(foamMesh->nCells(), 0.0);
    if (zoneI > 0) {
      for (int i = 0; i < zoneI; i++) {
        const Foam::cellZone &cz = cellZones[i];
        forAll (cz, j)
          physGrps[cz[j]] = i;
      }

      vtkSmartPointer<vtkDoubleArray> pg =
          vtkSmartPointer<vtkDoubleArray>::New();
      pg->SetName(phyGrpArrayName.c_str());
      pg->SetNumberOfComponents(1);
      for (int i = 0; i < foamMesh->nCells(); i++)
        pg->InsertNextTuple1(physGrps[i]);
      vtkdataSet->GetCellData()->AddArray(pg);
    }

    // Create boundary patches arrays and add to VTK
    const Foam::polyBoundaryMesh &patches = foamMesh->boundaryMesh();
    Foam::wordList patchNames = patches.names();
    Foam::wordList patchTypes = patches.types();
    // const Foam::labelList &patchIds = patches.patchID();
    Foam::label patchSize = patches.size();
    std::vector<int> face_patch_map;  // Which patch the current face is in
    Foam::label nFaces = foamMesh->nFaces();
    for (int i = 0; i < nFaces; i++)
      face_patch_map.push_back(int(patches.whichPatch(i)));

    // Create sideSets
    vtkSmartPointer<vtkPolyData> sideSet = vtkSmartPointer<vtkPolyData>::New();
    Foam::faceList allFaces = foamMesh->faces();
    Foam::labelList faceOwners = foamMesh->faceOwner();
    Foam::pointField allPoints = foamMesh->points();

    // Filter internal faces out
    int nonInternalFaces = 0;
    for (int i = 0; i < patchSize; i++) nonInternalFaces += patches[i].size();
    Foam::faceList sideSetOnlyFaces(nonInternalFaces);
    Foam::labelList sideSetOnlyFaceLabels(nonInternalFaces);
    Foam::labelList sideSetOnlyOwners(nonInternalFaces);
    int indx = 0;
    for (int i = 0; i < static_cast<int>(face_patch_map.size()); i++) {
      if (face_patch_map[i] != -1) {
        sideSetOnlyFaces[indx] = allFaces[i];
        sideSetOnlyOwners[indx] = faceOwners[i];
        sideSetOnlyFaceLabels[indx] = i;
        indx++;
      }
    }

    sideSet->SetPoints(vtkdataSet->GetPoints());
    sideSet->Allocate();

    // Start adding cells (2D faces)
    std::vector<int> facePntIds;
    for (int i = 0; i < sideSetOnlyFaces.size(); i++) {
      if (sideSetOnlyFaces[i].size() == 3) {
        // add triangle face (#5)
        vtkSmartPointer<vtkIdList> vtkCellIds =
            vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(sideSetOnlyFaces[i].size());
        for (int j = 0; j < sideSetOnlyFaces[i].size(); j++)
          vtkCellIds->SetId(j, sideSetOnlyFaces[i][j]);
        sideSet->InsertNextCell(5, vtkCellIds);
      } else if (sideSetOnlyFaces[i].size() == 4) {
        // add quad face (#9)
        vtkSmartPointer<vtkIdList> vtkCellIds =
            vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(sideSetOnlyFaces[i].size());
        for (int j = 0; j < sideSetOnlyFaces[i].size(); j++)
          vtkCellIds->SetId(j, sideSetOnlyFaces[i][j]);
        sideSet->InsertNextCell(9, vtkCellIds);
      } else if (sideSetOnlyFaces[i].size() > 4) {
        //  add polygon face (#7)
        vtkSmartPointer<vtkIdList> vtkCellIds =
            vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(sideSetOnlyFaces[i].size());
        for (int j = 0; j < sideSetOnlyFaces[i].size(); j++)
          vtkCellIds->SetId(j, sideSetOnlyFaces[i][j]);
        sideSet->InsertNextCell(7, vtkCellIds);
      } else {
        Foam::Info << "Face type not supported yet!" << Foam::nl;
      }
    }

    auto sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
    sideSetEntities->SetName("GeoEnt");

    auto sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
    sideSetOrigCellId->SetName("OrigCellIds");
    sideSetOrigCellId->SetNumberOfComponents(2);

    auto sideSetCellFaceId = vtkSmartPointer<vtkIntArray>::New();
    sideSetCellFaceId->SetName("CellFaceIds");
    sideSetCellFaceId->SetNumberOfComponents(2);

    auto sideSetPatchId = vtkSmartPointer<vtkIntArray>::New();
    sideSetPatchId->SetName("PatchIds");

    auto sideSetPatchName = vtkSmartPointer<vtkStringArray>::New();
    sideSetPatchName->SetName("PatchNames");

    for (int i = 0; i < sideSetOnlyFaces.size(); i++) {
      sideSetEntities->InsertNextValue(
          static_cast<int>(physGrps[sideSetOnlyOwners[i]]));
      sideSetOrigCellId->InsertNextTuple2(sideSetOnlyOwners[i],-1);
    }

    for (int i = 0; i < sideSetOnlyOwners.size(); i++) {
      const Foam::cell& facesCell = foamMesh->cells()[sideSetOnlyOwners[i]];
      for (int j = 0; j < facesCell.size(); j++) {
        if (foamMesh->faces()[facesCell[j]] == sideSetOnlyFaces[i]) {
          sideSetCellFaceId->InsertNextTuple2(j,-1);
        }
      }
    }

    for (int i = 0; i < sideSetOnlyFaces.size(); i++) {
      sideSetPatchId->InsertNextValue(face_patch_map[sideSetOnlyFaceLabels[i]]);
      sideSetPatchName->InsertNextValue(
          patchNames[face_patch_map[sideSetOnlyFaceLabels[i]]]);
    }

    // Add arrays to sideSet
    sideSet->GetCellData()->AddArray(sideSetPatchId);
    sideSet->GetFieldData()->AddArray(sideSetPatchName);

    // Create sideSet struct
    auto sideSetStruct =
        SideSet(sideSet, sideSetEntities, sideSetOrigCellId, sideSetCellFaceId);

    std::string gmshMesh = "foamGeoMesh_" + nemAux::getRandomString(6);
#ifdef HAVE_GMSH
    GmshInterface::Initialize();
    gmsh::model::add(gmshMesh);
    gmsh::model::setCurrent(gmshMesh);

    {  // Add geometric entities and physical groups
      std::set<std::pair<int, int>> dim_phyGrp;
      int idx = 0;
      auto phyGrpArray = vtkArrayDownCast<vtkDataArray>(
          vtkdataSet->GetCellData()->GetArray(phyGrpArrayName.c_str(), idx));
      vtkSmartPointer<vtkGenericCell> vtkGC =
          vtkSmartPointer<vtkGenericCell>::New();

      if (idx != -1) {
        // First sort
        for (vtkIdType i = 0; i < vtkdataSet->GetNumberOfCells(); ++i) {
          vtkGC->SetCellType(vtkdataSet->GetCellType(i));
          int dim = vtkGC->GetCellDimension();
          int phyGrp = static_cast<int>(phyGrpArray->GetComponent(i, 0));

          dim_phyGrp.insert({dim, phyGrp});
        }

        // then add. Each phyGrp gets its own geoEnt
        for (const auto &dp : dim_phyGrp) {
          gmsh::model::addDiscreteEntity(dp.first, dp.second);
          gmsh::model::addPhysicalGroup(dp.first, {dp.second}, dp.second);
        }
      }
    }
#endif
    return {vtkdataSet, gmshMesh, phyGrpArrayName, sideSetStruct};
  }
}

std::unique_ptr<Foam::fvMesh> foamGeoMesh::GM2foam(
    const GeoMesh &geoMesh, Foam::Time *runTime,
    const std::string &phyGrpArrayName) {
  int numPoints = static_cast<int>(geoMesh.mesh->GetNumberOfPoints());
  int numCells = static_cast<int>(geoMesh.mesh->GetNumberOfCells());

  Foam::pointField pointData(numPoints);
  Foam::pointField pointData2(numPoints);
  Foam::cellShapeList cellShapeData(numCells);
  Foam::cellShapeList cellShapeData2(numCells);

  if (numPoints > 0) {
    // Fetch all points and coordinates
    std::vector<std::vector<double>> verts;
    verts.resize(numPoints);
    for (int ipt = 0; ipt < numPoints; ipt++) {
      std::vector<double> getPt = std::vector<double>(3);
      geoMesh.mesh->GetPoint(ipt, &getPt[0]);
      verts[ipt].resize(3);
      verts[ipt][0] = getPt[0];
      verts[ipt][1] = getPt[1];
      verts[ipt][2] = getPt[2];
    }

    // Gets Ids for cells
    std::vector<std::vector<int>> cellIds;
    cellIds.resize(numCells);

    // Gets celltypes for all cells in mesh
    std::vector<int> typeCell;
    typeCell.resize(numCells);

    for (int i = 0; i < numPoints; i++) {
      pointData[i] = Foam::vector(verts[i][0], verts[i][1], verts[i][2]);
      pointData2[i] = Foam::vector(verts[i][0], verts[i][1], verts[i][2]);
    }

    for (int i = 0; i < numCells; i++) {
      vtkIdList *ptIds = vtkIdList::New();
      geoMesh.mesh->GetCellPoints(i, ptIds);
      int numIds = static_cast<int>(ptIds->GetNumberOfIds());
      cellIds[i].resize(numIds);
      for (int j = 0; j < numIds; j++) cellIds[i][j] = static_cast<int>(ptIds->GetId(j));
      Foam::labelList meshPoints(numIds);
      for (int k = 0; k < numIds; k++) meshPoints[k] = cellIds[i][k];
      typeCell[i] = geoMesh.mesh->GetCellType(i);
      if (typeCell[i] == 12) {
        cellShapeData[i] = Foam::cellShape("hex", meshPoints, true);
        cellShapeData2[i] = Foam::cellShape("hex", meshPoints, true);
      } else if (typeCell[i] == 14) {
        cellShapeData[i] = Foam::cellShape("pyr", meshPoints, true);
        cellShapeData2[i] = Foam::cellShape("pyr", meshPoints, true);
      } else if (typeCell[i] == 10) {
        cellShapeData[i] = Foam::cellShape("tet", meshPoints, true);
        cellShapeData2[i] = Foam::cellShape("tet", meshPoints, true);
      } else {
        std::cerr << "Only Hexahedral, Tetrahedral,"
                  << " and Pyramid cells are supported in foamGeoMesh!"
                  << std::endl;
        throw;
      }
    }
  } else {
    return nullptr;
  }

  // Foam::word regName = Foam::polyMesh::defaultRegion;
  Foam::word regName = "";

  // Get patch names from mesh
  int idx1 = 0;
  int idx2 = 0;
  std::map<int, std::string> ptchNameIDMap;
  Foam::faceListList bndryFaces(0);
  Foam::wordList BndryPatchNames(0);
  Foam::wordList BndryPatchTypes(0);
  Foam::wordList boundaryPatchPhysicalTypes(0);
  int totalPatches;
  std::vector<int> startIds;
  std::vector<int> nFacesInPatch;
  Foam::PtrList<Foam::polyPatch> patchPtrList(0);
  Foam::cellList cellLst2(0);
  Foam::faceList faceLst2(0);

  if (geoMesh.sideSet.sides) {
    auto ptchIds = vtkIntArray::FastDownCast(
        geoMesh.sideSet.sides->GetCellData()->GetArray("PatchIds", idx1));
    auto ptchNms = vtkStringArray::SafeDownCast(
        geoMesh.sideSet.sides->GetFieldData()->GetAbstractArray("PatchNames",
                                                                idx2));

    if (idx1 != -1 && idx2 != -1) {
      // Patches exist
      for (int i = 0; i < ptchIds->GetNumberOfTuples(); ++i) {
        // Get patch Id
        auto cc = ptchIds->GetTypedComponent(i, 0);
        int id_p = cc;

        // Get patch name
        std::string nm = ptchNms->GetValue(i);

        // Create maps
        ptchNameIDMap[id_p] = nm;
      }

      totalPatches = static_cast<int>(ptchNameIDMap.size());
      std::vector<int> faceCounts(totalPatches,
                                  0);  // Counts faces for each patch
      std::vector<std::vector<int>> facePatchIds;
      facePatchIds.resize(totalPatches);
      for (int i = 0; i < ptchIds->GetNumberOfTuples(); ++i) {
        // Get patch Id
        auto cc = ptchIds->GetTypedComponent(i, 0);
        int id_p = cc;
        faceCounts[id_p] += 1;
        facePatchIds[id_p].push_back(i);
      }

      BndryPatchNames.setSize(totalPatches);
      BndryPatchTypes.setSize(totalPatches);

      auto iter = ptchNameIDMap.begin();
      int r_indx = 0;
      while (iter != ptchNameIDMap.end()) {
        BndryPatchNames[r_indx] = iter->second;
        BndryPatchTypes[r_indx] = "patch";
        iter++;
        r_indx++;
      }

      // Start pulling out cell (2D faces) point order from sideSets and create
      // faceList. Append those into faceListList at the end of every patch.
      bndryFaces.clear();
      for (int i = 0; i < totalPatches; i++) {
        Foam::faceList thisPatch(0);
        thisPatch.clear();
        for (int j : facePatchIds[i]) {
          vtkIdList *ptIds = vtkIdList::New();
          geoMesh.sideSet.sides->GetCellPoints(j, ptIds);
          int numIds = static_cast<int>(ptIds->GetNumberOfIds());
          Foam::labelList lstLbl(0);
          lstLbl.clear();
          for (int k = 0; k < numIds; k++) {
            lstLbl.append(static_cast<int>(ptIds->GetId(k)));
          }
          if (!lstLbl.empty()) thisPatch.append(Foam::face(lstLbl));
        }
        bndryFaces.append(thisPatch);
        nFacesInPatch.push_back(thisPatch.size());
      }
      boundaryPatchPhysicalTypes.resize(totalPatches,
                                        Foam::polyPatch::typeName);
    }
  }

  {
    Foam::polyMesh tmpfm(
        Foam::IOobject(regName, runTime->timeName(), *runTime,
                       Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
        std::move(pointData),       // Vertices
        cellShapeData,              // Cell shape and points
        bndryFaces,                 // Boundary faces
        BndryPatchNames,            // Boundary Patch Names
        BndryPatchTypes,            // Boundary Patch Types
        "defaultPatch",             // Default Patch Name
        Foam::polyPatch::typeName,  // Default Patch Type
        Foam::wordList());

    const Foam::polyBoundaryMesh &patches2 = tmpfm.boundaryMesh();

    for (int i = 0; i < patches2.size(); i++) {
      auto *tmpPtch = new Foam::polyPatch(patches2[i]);
      patchPtrList.append(tmpPtch);
    }

    // Get pointField (points), faceList (faces), and cellList (cells).
    const Foam::cellList &cellLst = tmpfm.cells();
    cellLst2 = cellLst;
    const Foam::faceList &faceLst = tmpfm.faces();
    faceLst2 = faceLst;
  }

  auto fm = std::unique_ptr<Foam::fvMesh>(new Foam::fvMesh(
      Foam::IOobject(regName, runTime->timeName(), *runTime,
                     Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
      std::move(pointData2), std::move(faceLst2), std::move(cellLst2)));

  // Add boundary patches to fvMesh using patchPtrList
  fm->addFvPatches(patchPtrList, true);

  // Get Physical Groups Vector from GeoMesh
  std::vector<double> allGroups;
  int idx;
  vtkSmartPointer<vtkDataArray> cd =
      geoMesh.mesh->GetCellData()->GetArray(phyGrpArrayName.c_str(), idx);
  if (idx != -1) {
    // allGroups.resize(cd->GetNumberOfTuples());
    double x[1];
    for (nemId_t i = 0; i < cd->GetNumberOfTuples(); ++i) {
      cd->GetTuple(i, x);
      allGroups.push_back(x[0]);
    }

    // Figure out number of physical groups
    std::vector<double> scratchVec = allGroups;
    std::sort(scratchVec.begin(), scratchVec.end());
    scratchVec.erase(unique(scratchVec.begin(), scratchVec.end()),
                     scratchVec.end());
    int totalGroups = static_cast<int>(scratchVec.size());

    // Create cellZones
    Foam::label zoneI;
    for (int i = 0; i < totalGroups; i++) {
      Foam::labelList currentGroupList;
      //for (int j : allGroups) {
      for (int j = 0; j < static_cast<int>(allGroups.size()); j++) {
        if (static_cast<int>(allGroups[j]) == i) currentGroupList.append(j);
      }

      const Foam::cellZoneMesh &cellZones = fm->cellZones();
      zoneI = cellZones.size();
      fm->cellZones().setSize(zoneI + 1);
      Foam::word nameSet = "pg_" + std::to_string(i);
      auto clzn =
          new Foam::cellZone(nameSet, currentGroupList, zoneI, cellZones);
      fm->cellZones().set(zoneI, clzn);
    }
    fm->cellZones().writeOpt() = Foam::IOobject::AUTO_WRITE;
  }
  return fm;
}

void foamGeoMesh::write(const std::string &fileName) {
  // Create system directory
  Foam::word w = "system";
  if (!Foam::isDir(w)) Foam::mkDir(w);

  // Write ControlDict
  Foam::fileName fcontrolDict = "system/controlDict";
  if (!Foam::exists(fcontrolDict)) {
    Foam::OFstream outcontrolDict(fcontrolDict);
    Foam::IOobject::writeBanner(outcontrolDict);
    controlDict_->write(outcontrolDict, false);
  }

  // Write fvSchemes
  Foam::fileName ffvSchemes = "system/fvSchemes";
  if (!Foam::exists(ffvSchemes)) {
    Foam::OFstream outfvSchemes(ffvSchemes);
    Foam::IOobject::writeBanner(outfvSchemes);
    fvSchemes_->write(outfvSchemes, false);
  }

  // Write fvSolution
  Foam::fileName ffvSolution = "system/fvSolution";
  if (!Foam::exists(ffvSolution)) {
    Foam::OFstream outfvSolution(ffvSolution);
    Foam::IOobject::writeBanner(outfvSolution);
    fvSolution_->write(outfvSolution, false);
  }

  // Region to write in
  Foam::fileName newReg = "";
  if (fileName.empty())
    newReg = "";
  else
    newReg = fileName;

  // Write dictionaries to the specific region too
  if (!Foam::isDir("system/" + newReg)) Foam::mkDir("system/" + newReg);

  Foam::fileName fcontrolDict_2 = "system/" + newReg + "/controlDict";
  if (!Foam::exists(fcontrolDict_2)) { Foam::cp(fcontrolDict, fcontrolDict_2); }

  // Write fvSchemes
  Foam::fileName ffvSchemes_2 = "system/" + newReg + "/fvSchemes";
  if (!Foam::exists(ffvSchemes_2)) { Foam::cp(ffvSchemes, ffvSchemes_2); }

  // Write fvSolution
  Foam::fileName ffvSolution_2 = "system/" + newReg + "/fvSolution";
  if (!Foam::exists(ffvSolution_2)) { Foam::cp(ffvSolution, ffvSolution_2); }

  // Write mesh in fileName
  if (!Foam::isDir("0")) Foam::mkDir("0");

  // Write cellZones and cellSets
  const Foam::cellZoneMesh &czMesh = fmesh_->cellZones();
  Foam::label numZones = czMesh.size();

  fmesh_->setInstance("0");
  fmesh_->write();

  if (numZones > 0) {
    for (int i = 0; i < numZones; i++) {
      const Foam::cellZone &cz = czMesh[i];
      Foam::cellSet(*fmesh_, cz.name(), cz).write();
    }
  }

  Foam::fileName currentDir = "0";
  Foam::fileName currentReg = fmesh_->dbDir();
  Foam::fileName newDir = "constant";
  Foam::fileName finalMesh = "";

  if (currentReg.empty()) {
    if (!Foam::isDir(currentDir + "/" + newReg))
      Foam::mkDir(currentDir + "/" + newReg);
    if (!newReg.empty()) {
      Foam::mv(currentDir + "/polyMesh",
               currentDir + "/" + newReg + "/polyMesh");
      finalMesh = currentDir + "/" + newReg;
    } else {
      finalMesh = currentDir + "/polyMesh";
    }
  } else {
    if (!newReg.empty()) {
      Foam::mv(currentDir + "/" + currentReg, currentDir + "/" + newReg);
      finalMesh = currentDir + "/" + newReg;
    } else {
      Foam::mv(currentDir + "/" + currentReg, currentDir + "/");
      finalMesh = currentDir + "/polyMesh";
    }
  }

  if (!Foam::isDir(newDir + "/" + newReg)) {
    if (!newReg.empty()) {
      Foam::mkDir(newDir + "/" + newReg);
      Foam::mv(finalMesh, newDir + "/" + newReg);
    } else {
      Foam::mkDir(newDir + "/polyMesh");
      Foam::mv(finalMesh, newDir + "/polyMesh");
    }
  }

  if (!newReg.empty()) {
    if (!Foam::isDir(newDir + "/" + newReg)) {
      Foam::mkDir(newDir + "/" + newReg);
      Foam::mv(finalMesh, newDir + "/" + newReg);
    }
  } else {
    if (!Foam::isDir(newDir + "/polyMesh")) Foam::mkDir(newDir + "/polyMesh");
    Foam::mv(finalMesh, newDir + "/polyMesh");
  }
}

void foamGeoMesh::report(std::ostream &out) const { geoMeshBase::report(out); }

void foamGeoMesh::setFoamMesh(std::unique_ptr<Foam::fvMesh> foamMesh) {
  // Setting class parameters
  fmesh_ = std::unique_ptr<Foam::fvMesh>(foamMesh.get());

  // Parent
  setGeoMesh(foam2GM(foamMesh.get()));
}

void foamGeoMesh::resetNative() {
  // Clear
  controlDict_.reset();
  fvSchemes_.reset();
  fvSolution_.reset();
  runTime_.reset();
  fmesh_.reset();

  // Initialize required variables
  InitializeFoam();

  // Reset fvMesh
  fmesh_ = GM2foam(getGeoMesh(), runTime_.get());
}

void foamGeoMesh::InitializeFoam() {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  controlDict_ = initFoam->createControlDict(writeDicts);
  fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  runTime_ =
      std::unique_ptr<Foam::Time>(new Foam::Time(controlDict_.get(), ".", "."));
  Foam::argList::noParallel();
}

}  // namespace MSH
}  // namespace NEM

