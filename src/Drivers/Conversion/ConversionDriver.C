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
#include "Drivers/Conversion/ConversionDriver.H"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <utility>

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkModelMetadata.h>
#include <vtkPolyData.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#include "Mesh/cobalt.H"
#include "MeshOperation/meshSrch.H"
#include "Mesh/geoMeshFactory.H"

namespace NEM {
namespace DRV {

void ConversionDriver::genExo(meshBase *mb, NEM::MSH::EXOMesh::exoMesh *em,
                              const int &ndeIdOffset, const int &elmIdOffset,
                              int &ins, int &ieb, int &iss, std::string mshName,
                              const bool &usePhys, int &ndeIdOffset_local,
                              int &elmIdOffset_local,
                              const bool &makeFreeSurfSS,
                              const bool &splitTopBotSS,
                              std::vector<std::string> sideSetNames) {
  // add nodes to database
  for (nemId_t iNde = 0; iNde < mb->getNumberOfPoints(); ++iNde)
    em->addNde(mb->getPoint(iNde));

  // node coordinate to one nodeSet
  NEM::MSH::EXOMesh::ndeSetType ns;
  ns.id = ++ins;
  ns.nNde = mb->getNumberOfPoints();
  ns.name = mshName;
  ns.ndeIdOffset = ndeIdOffset;
  for (int iNde = 0; iNde < ns.nNde; iNde++) {
    ns.ndeIds.emplace_back(iNde + 1);
  }
  em->addNdeSet(ns);
  ndeIdOffset_local += ns.nNde;

  // add element blocks

  // Element bucket where they will be sorted.
  // Outer layer: Specifies grouping, such as physical groups. The 0th index
  // is reserved for elements without a group.
  // Inner layer: Sorted by element type. The 0th index is reserved for
  // unsupported elements. The current support respects the ordering of the
  // NEM::MSH::EXOMesh::elementType enum:
  //     OTHER, TRIANGLE, QUAD, TETRA, HEX
  std::map<int, std::map<NEM::MSH::EXOMesh::elementType, std::vector<int>>>
      elmBucket;
  // map for VTK to EXO element ids
  std::map<int, int> v2e_elemID_map;

  std::vector<double> grpIds(mb->getNumberOfCells(), 0.0);
  if (usePhys) mb->getCellDataArray("PhysGrpId", grpIds);

  bool is3D = false;
  // Loop through all elements
  for (int iElm = 0; iElm < mb->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mb->getDataSet()->GetCellType(iElm));
    NEM::MSH::EXOMesh::elementType exoType =
        NEM::MSH::EXOMesh::v2eEMap(vtkType);

    // set the dimension of Exo database based on present elements
    if (!is3D) {
      if (exoType == NEM::MSH::EXOMesh::elementType::TETRA ||
          exoType == NEM::MSH::EXOMesh::elementType::WEDGE ||
          exoType == NEM::MSH::EXOMesh::elementType::HEX) {
        is3D = true;
        em->setDimension(3);
      }
    }

    elmBucket[static_cast<int>(grpIds[iElm])][exoType].emplace_back(iElm);
  }

  // sanity check
  int numUnsupported = 0;
  for (const auto &elmGroup : elmBucket)
    if (elmGroup.second.count(NEM::MSH::EXOMesh::elementType::OTHER) != 0)
      numUnsupported +=
          elmGroup.second.at(NEM::MSH::EXOMesh::elementType::OTHER).size();
  if (numUnsupported > 0) {
    std::cerr << "WARNING: Detected " << numUnsupported
              << " unsupported elements.\n";
    throw;
  }

  // for each group and supported type, if existent, add an element block
  for (const auto &elmGroup : elmBucket) {
    for (const auto &elmIds : elmGroup.second) {
      if (elmIds.second.empty()) continue;  // skip if empty

      NEM::MSH::EXOMesh::elmBlkType eb;
      eb.id = ++ieb;
      eb.ndeIdOffset = ndeIdOffset;
      eb.nElm = elmIds.second.size();

      // populate VTK to EXO element id map for side set implementation
      for (int i = 0; i < elmIds.second.size(); ++i) {
        int vtkId = elmIds.second[i];
        int exoId = i + elmIdOffset_local + 1;
        v2e_elemID_map.insert(std::pair<int, int>(vtkId, exoId));
      }

      if (usePhys)
        eb.name = mshName + "_PhysGrp_" + std::to_string(ieb);
      else
        eb.name = mshName + "_" + std::to_string(ieb);

      switch (elmIds.first) {
        case NEM::MSH::EXOMesh::elementType::TRIANGLE:
          std::cout << "Number of triangular elements = " << eb.nElm << "\n";
          eb.ndePerElm = 3;
          eb.eTpe = NEM::MSH::EXOMesh::elementType::TRIANGLE;
          break;
        case NEM::MSH::EXOMesh::elementType::QUAD:
          std::cout << "Number of quadrilateral elements = " << eb.nElm << "\n";
          eb.ndePerElm = 4;
          eb.eTpe = NEM::MSH::EXOMesh::elementType::QUAD;
          break;
        case NEM::MSH::EXOMesh::elementType::TETRA:
          std::cout << "Number of tetrahedral elements = " << eb.nElm << "\n";
          eb.ndePerElm = 4;
          eb.eTpe = NEM::MSH::EXOMesh::elementType::TETRA;
          break;
        case NEM::MSH::EXOMesh::elementType::HEX:
          std::cout << "Number of hexahedral elements = " << eb.nElm << "\n";
          eb.ndePerElm = 8;
          eb.eTpe = NEM::MSH::EXOMesh::elementType::HEX;
          break;
        case NEM::MSH::EXOMesh::elementType::WEDGE:
          std::cout << "Number of wedge elements = " << eb.nElm << "\n";
          eb.ndePerElm = 6;
          eb.eTpe = NEM::MSH::EXOMesh::elementType::WEDGE;
          break;
        case NEM::MSH::EXOMesh::elementType::OTHER:
        default:
          std::cerr << "WARNING: Processing unsupported element. Previous "
                       "sanity check failed!\n";
          throw;
      }

      eb.conn.reserve(eb.nElm * eb.ndePerElm);
      vtkIdList *nids = vtkIdList::New();
      for (const auto &iElm : elmIds.second) {
        mb->getDataSet()->GetCellPoints(iElm, nids);
        for (int in = 0; in < eb.ndePerElm; ++in) {
          // offset node ids by 1
          eb.conn.emplace_back(nids->GetId(in) + 1);
        }
      }

      // DEBUG
      // std::cout << "Min node index = "
      //           << *min_element(eb.conn.begin(), eb.conn.end()) << "\n"
      //           << "Max node index = "
      //           << *max_element(eb.conn.begin(), eb.conn.end()) << "\n";
      // std::cout << "Starting node offset = " << ndeIdOffset << std::endl;

      em->addElmBlk(eb);
      elmIdOffset_local += eb.nElm;
    }
  }

  // add side set of the external surface(s) (free surface)
  if (makeFreeSurfSS) {
    freeSurfaceSideSet(mb, em, elmIdOffset, v2e_elemID_map, splitTopBotSS,
                       sideSetNames);
  }

  // add side sets
  if (mb->getMetadata()) {
    // get side set metadata
    vtkSmartPointer<vtkModelMetadata> metadata = mb->getMetadata();
    vtkSmartPointer<vtkStringArray> sdeSetNames = metadata->GetSideSetNames();
    int *sdeSetElmLst = metadata->GetSideSetElementList();
    int *sdeSetSdeLst = metadata->GetSideSetSideList();
    int *sdeSetSze = metadata->GetSideSetSize();

    for (int iSS = 0; iSS < metadata->GetNumberOfSideSets(); iSS++) {
      NEM::MSH::EXOMesh::sdeSetType ss;
      ss.id = ++iss;
      ss.name = sdeSetNames->GetValue(iSS);
      ss.nSde = sdeSetSze[iSS];
      ss.elmIds.assign(sdeSetElmLst, sdeSetElmLst + sdeSetSze[iSS]);
      ss.sdeIds.assign(sdeSetSdeLst, sdeSetSdeLst + sdeSetSze[iSS]);
      ss.elmIdOffset = elmIdOffset;
      em->addSdeSet(ss);

      // Advance pointer for reading next side set.
      sdeSetElmLst += sdeSetSze[iSS];
      sdeSetSdeLst += sdeSetSze[iSS];
    }
  }
}
void ConversionDriver::genExo(std::vector<meshBase *> meshes,
                              const std::string &fname) {
  bool usePhys = false;
  bool makeFreeSurfSS = false;
  bool splitTopBotSS = false;
  std::vector<std::string> sideSetNames;

  std::cout << "Warning: this method does not support physical groups."
            << std::endl;
  std::cout << "Warning: this method does not support post-processing."
            << std::endl;

  int nMsh = meshes.size();

  // sanity check
  if (nMsh == 0) {
    std::cerr << "Error: At least one mesh should be provided!\n";
    exit(-1);
  }

  // starting conversion operation
  auto em = new NEM::MSH::EXOMesh::exoMesh(fname);

  // reading meshes
  int ndeIdOffset = 0;
  int elmIdOffset = 0;
  int ins = 0;
  int ieb = 0;
  int iss = 0;
  int ndeIdOffset_local = 0;
  int elmIdOffset_local = 0;

  std::string mshName;

  for (auto itrMsh = meshes.begin(); itrMsh != meshes.end(); itrMsh++) {
    mshName = (*itrMsh)->getFileName();
    genExo(*itrMsh, em, ndeIdOffset, elmIdOffset, ins, ieb, iss, mshName,
           usePhys, ndeIdOffset_local, elmIdOffset_local, makeFreeSurfSS,
           splitTopBotSS, sideSetNames);
  }

  // writing the file
  em->write();
  em->report();

  em->mergeNodes();

  // clean up
  delete em;
}

void ConversionDriver::freeSurfaceSideSet(
    const meshBase *mb, NEM::MSH::EXOMesh::exoMesh *em, int elmIdOffset,
    std::map<int, int> v2e_elemID_map, bool splitTopBotSS,
    std::vector<std::string> sideSetNames) {
  std::cout << "Creating side set from free surface." << std::endl;

  std::vector<std::pair<int, int>> freeSurfCellEdge;
  std::vector<std::pair<int, int>> hexFreeSurfCellFace;
  std::vector<std::pair<int, int>> wedgeFreeSurfCellFace;
  std::vector<std::pair<int, int>> tetraFreeSurfCellFace;
  std::vector<int> elementIds1, elementIds2, elementIds3, sideIds1, sideIds2,
      sideIds3;
  std::vector<std::vector<int>> splitElemIds, splitSideIds;

  // Maps for VTK to EXO face id conversion for hex and wedge
  std::map<int, int> v2e_hexFace_map = {{1, 4}, {2, 2}, {3, 1},
                                        {4, 3}, {5, 5}, {6, 6}};
  std::map<int, int> v2e_wedgeFace_map = {
      {1, 4}, {2, 5}, {3, 1}, {4, 2}, {5, 3}};
  std::map<int, int> v2e_tetraFace_map = {{1, 1}, {2, 2}, {3, 3}, {4, 4}};

  // acquiring dataset
  vtkSmartPointer<vtkDataSet> ds = mb->getDataSet();
  if (!ds) {
    std::cerr << "No dataset is associated to the meshbase." << std::endl;
    throw;
  }
  bool tetWarning = false;
  // loop through cells and obtain different quantities needed
  int nCl = ds->GetNumberOfCells();
  for (int cell_i = 0; cell_i < nCl; cell_i++) {
    vtkCell *vc = ds->GetCell(cell_i);
    if (vc->GetCellDimension() == 2) {
      // for 2D cells
      int ne = vc->GetNumberOfEdges();
      for (int edge_i = 0; edge_i < ne; edge_i++) {
        vtkCell *ve = vc->GetEdge(edge_i);
        vtkIdList *pidl = ve->GetPointIds();

        std::pair<int, int> adjPair;

        // getting neighbour list
        vtkSmartPointer<vtkIdList> cidl = vtkSmartPointer<vtkIdList>::New();
        ds->GetCellNeighbors(cell_i, pidl, cidl);
        if (cidl->GetNumberOfIds() == 0) {
          adjPair.first = cell_i;
          adjPair.second = edge_i + 1;
          freeSurfCellEdge.push_back(adjPair);
        }
      }
    } else if (vc->GetCellDimension() == 3) {
      // for 3D cells
      int nfc = vc->GetNumberOfFaces();
      for (int face_i = 0; face_i < nfc; face_i++) {
        vtkCell *vf = vc->GetFace(face_i);
        vtkIdList *facePoints = vf->GetPointIds();

        std::pair<int, int> adjPair;

        // getting neighbour list
        vtkSmartPointer<vtkIdList> cidl = vtkSmartPointer<vtkIdList>::New();
        ds->GetCellNeighbors(cell_i, facePoints, cidl);
        if (cidl->GetNumberOfIds() == 0) {
          adjPair.first = cell_i;
          adjPair.second = face_i + 1;
          if (nfc == 6) hexFreeSurfCellFace.push_back(adjPair);
          if (nfc == 5) wedgeFreeSurfCellFace.push_back(adjPair);
          if (nfc == 4) {
            tetraFreeSurfCellFace.push_back(adjPair);
            if (splitTopBotSS && tetWarning == false) {
              std::cerr
                  << "Warning: Tetrahedral element found on free surface.\n"
                  << "         Cannot split top and bottom into separate side "
                     "sets.\n"
                  << "         Creating single side set." << std::endl;
              splitTopBotSS = false;
              tetWarning = true;
            }
          }
        }
      }
    }
  }

  if (em->getDimension() == 3) {
    // Iterate through *freeSurfCellFace to make lists of elements/faces
    // Hexahedra
    for (int i = 0; i < hexFreeSurfCellFace.size(); i++) {
      int exoId = v2e_elemID_map[hexFreeSurfCellFace[i].first];
      int sId = v2e_hexFace_map[hexFreeSurfCellFace[i].second];

      if (splitTopBotSS) {
        if (sId == 5) {
          // Top
          elementIds2.push_back(exoId);
          sideIds2.push_back(sId);
        } else if (sId == 6) {
          // Bottom
          elementIds3.push_back(exoId);
          sideIds3.push_back(sId);
        } else {
          elementIds1.push_back(exoId);
          sideIds1.push_back(sId);
        }
      } else {
        elementIds1.push_back(exoId);
        sideIds1.push_back(sId);
      }
    }
    // Wedges
    for (int i = 0; i < wedgeFreeSurfCellFace.size(); i++) {
      int exoId = v2e_elemID_map[wedgeFreeSurfCellFace[i].first];
      int sId = v2e_wedgeFace_map[wedgeFreeSurfCellFace[i].second];

      if (splitTopBotSS) {
        if (sId == 4) {
          // Top
          elementIds2.push_back(exoId);
          sideIds2.push_back(sId);
        } else if (sId == 5) {
          // Bottom
          elementIds3.push_back(exoId);
          sideIds3.push_back(sId);
        } else {
          elementIds1.push_back(exoId);
          sideIds1.push_back(sId);
        }
      } else {
        elementIds1.push_back(exoId);
        sideIds1.push_back(sId);
      }
    }
    // Tetrahedra
    for (int i = 0; i < tetraFreeSurfCellFace.size(); i++) {
      int exoId = v2e_elemID_map[tetraFreeSurfCellFace[i].first];
      int sId = v2e_tetraFace_map[tetraFreeSurfCellFace[i].second];
      elementIds1.push_back(exoId);
      sideIds1.push_back(sId);
    }
  }

  if (em->getDimension() == 2) {
    // Iterate through freeSurfCellEdge to make lists of elements/faces
    for (int i = 0; i < freeSurfCellEdge.size(); i++) {
      // map from VTK element id to EXO element id
      int exoId = v2e_elemID_map[freeSurfCellEdge[i].first];
      elementIds1.push_back(exoId);
      sideIds1.push_back(freeSurfCellEdge[i].second);
    }
  }

  if (splitTopBotSS && em->getDimension() == 3) {
    splitElemIds = {elementIds1, elementIds2, elementIds3};
    splitSideIds = {sideIds1, sideIds2, sideIds3};
  } else {
    splitElemIds = {elementIds1};
    splitSideIds = {sideIds1};
  }

  for (int i = 0; i < splitElemIds.size(); ++i) {
    NEM::MSH::EXOMesh::sdeSetType ss;

    // TODO: fix names, still not working correctly
    if (sideSetNames.size() == 1 && splitElemIds.size() == 1) {
      ss.name = sideSetNames[i];
    }
    if (sideSetNames.size() == 0) {
      std::cout << "side set names size is zero" << std::endl;
      ss.name = "side_set_000000" + std::to_string(i + 1);
    }
    if (sideSetNames.size() > 0 && sideSetNames.size() < 3 && splitTopBotSS) {
      std::cerr << "Error: Expected 3 'Side Set Names', found "
                << sideSetNames.size() << std::endl;
      exit(1);
    }
    if (sideSetNames.size() == 3) {
      ss.name = sideSetNames[i];
    }

    // ss.id = ++iss;
    ss.id = i + 1;
    ss.nSde = (int)splitSideIds[i].size();
    ss.elmIds = splitElemIds[i];
    ss.sdeIds = splitSideIds[i];
    ss.elmIdOffset = elmIdOffset;
    em->addSdeSet(ss);
  }
}

void ConversionDriver::procExo(const jsoncons::json &ppJson,
                               const std::string &fname,
                               NEM::MSH::EXOMesh::exoMesh *em) {
  // converting to mesh base for geometric inquiry
  meshBase *mb = meshBase::Create(fname);

  // performing requested operation
  std::string opr = ppJson.get_with_default("Operation", "");
  if (opr == "Material Assignment") {
    meshSrch *ms = meshSrch::Create(mb);
    // gathering information about all zones
    // if densities are defined, materials with higher density will be
    // prioritized
    bool appDen = ppJson.get_with_default("Apply Density", false);

    std::map<std::pair<double, std::string>, std::set<int>> zoneGeom;

    // loop over all zones
    for (const auto &zone : ppJson["Zones"].array_range()) {
      // assuming first element is zone information keyed by zone name
      // that we do not care about yet
      std::string matName = zone[0].get_with_default("Material Name", "N/A");
      std::string shape = zone[0].get_with_default("Shape", "N/A");
      std::cout << "Processing zone: " << zone.object_range().begin()->key()
                << "  Material: " << matName << "  Shape: " << shape;

      double density = 1.0;  // Default density is 1.0
      if (appDen) {
        density = zone[0].get_with_default("Density", 1.0);
        std::cout << "  Density: " << density;
      }
      std::cout << std::endl;

      std::vector<nemId_t> lst;
      if (shape == "Box") {
        // Box shape. Requires 3-vector of Min and Max in x, y, and z, resp.
        std::vector<double> bb;
        bb.push_back(zone[0]["Params"]["Min"][0].as<double>());
        bb.push_back(zone[0]["Params"]["Max"][0].as<double>());
        bb.push_back(zone[0]["Params"]["Min"][1].as<double>());
        bb.push_back(zone[0]["Params"]["Max"][1].as<double>());
        bb.push_back(zone[0]["Params"]["Min"][2].as<double>());
        bb.push_back(zone[0]["Params"]["Max"][2].as<double>());

        ms->FindCellsWithinBounds(bb, lst, true);
      } else if (shape == "STL") {
        // STL shape. Only supports Tri surface.
        // Node Coordinates are given as 3-vector in an array.
        // Connectivities are given as 3-vectors in an array.
        // Tris are 0-indexed.
        std::vector<std::vector<double>> crds;
        std::vector<std::vector<vtkIdType>> conns;
        for (const auto &crd :
             zone[0]["Params"]["Node Coordinates"].array_range())
          crds.push_back(crd.as<std::vector<double>>());
        for (const auto &conn :
             zone[0]["Params"]["Connectivities"].array_range())
          conns.push_back(conn.as<std::vector<vtkIdType>>());

        ms->FindCellsInTriSrf(crds, conns, lst);
      } else if (shape == "Sphere") {
        // Sphere shape.
        // Center is a 3-vector in an array.
        // Radius is a double.
        std::vector<double> center =
            zone[0]["Params"]["Center"].as<std::vector<double>>();
        auto radius = zone[0]["Params"]["Radius"].as<double>();

        ms->FindCellsInSphere(center, radius, lst);
      } else {
        std::cerr << "WARNING: Skipping unknown zone shape: " << shape
                  << std::endl;
        continue;
      }

      if (zone[0].contains("Only From Block")) {
        std::string blkName = zone[0]["Only From Block"].as<std::string>();

        std::vector<std::string> elmBlkNames = em->getElmBlkNames();
        auto elmBlkName =
            std::find(elmBlkNames.begin(), elmBlkNames.end(), blkName);
        if (elmBlkName == elmBlkNames.end()) {
          std::cerr << "WARNING: Only From Block " << blkName
                    << " matches no available blocks. Continuing with no "
                       "restriction.\n";
        } else {
          std::vector<int> lst_int(lst.begin(), lst.end());
          bool allIn = false;
          lst_int = em->lstElmInBlk(
              std::distance(elmBlkNames.begin(), elmBlkName), lst_int, allIn);
          lst.assign(lst_int.begin(), lst_int.end());
        }
      }
      zoneGeom[{1.0 / density, matName}].insert(lst.begin(), lst.end());
    }

    if (appDen)
      std::cout << "Applying material zones based on density ordering"
                << std::endl;

    // adjusting exodus database accordingly
    for (const auto &zone : zoneGeom) {
      // zone is ((density, material name), ids)
      std::vector<int> elmLst;
      std::cout << "Manipulating ExodusDB for " << zone.first.second
                << std::endl;
      elmLst.insert(elmLst.end(), zone.second.begin(), zone.second.end());
      em->addElmBlkByElmIdLst(zone.first.second, elmLst);
    }
  } else if (opr == "Check Duplicate Elements") {
    std::cout << "Checking for existence of duplicate elements ... ";
    meshSrch *ms = meshSrch::Create(mb);
    bool ret = ms->chkDuplElm();
    if (ret) {
      std::cerr << " The exodus database contains duplicate elements."
                << std::endl;
      exit(-1);
    } else {
      std::cout << "False" << std::endl;
    }
  } else if (opr == "Remove Block") {
    std::string blkName = ppJson.get_with_default("Block Name", "");
    std::cout << "Removing Block " << blkName << std::endl;
    em->removeElmBlkByName(blkName);
  } else if (opr == "Snap Node Coords To Zero") {
    double tol = ppJson.get_with_default("Tolerance", 0.0);
    std::cout << "Snapping nodal coordinates to zero using tolerance: " << tol
              << std::endl;
    em->snapNdeCrdsZero(tol);
  } else if (opr == "Boundary Condition Assignment") {
    // For EP16 boundary conditions are simply translated to node sets. Node
    // sets may have shared nodes. In that case a node the order of nodeset
    // matter. A later node set supersedes an earlier one.
    meshSrch *ms = meshSrch::Create(mb);

    // gathering information about all boundary node sets
    jsoncons::json bcs = ppJson["Condition"];
    for (const auto &bc : bcs.array_range()) {
      std::set<nemId_t> pntIds;

      // identify node ids on each boundary
      std::string bcName = bc["Name"].as<std::string>();
      std::string bcTyp = bc["Boundary Type"].as<std::string>();
      if (bcTyp == "Faces") {
        std::vector<double> srfCrd;
        std::vector<nemId_t> srfConn;
        jsoncons::json nc = bc["Params"]["Node Coordinates"];
        for (const auto &crds : nc.array_range())
          for (const auto &cmp : crds.array_range())
            srfCrd.push_back(cmp.as<double>());
        jsoncons::json conn = bc["Params"]["Connectivities"];
        for (const auto &tri : conn.array_range())
          for (const auto &idx : tri.array_range())
            srfConn.push_back(idx.as<nemId_t>());
        ms->FindPntsOnTriSrf(srfCrd, srfConn, pntIds);
        std::cout << "Number of points residing on the boundary " << bcName
                  << " is " << pntIds.size() << std::endl;
      } else if (bcTyp == "Edges") {
        std::vector<double> edgeCrd;
        jsoncons::json ncs = bc["Params"]["Start"];
        for (const auto &crds : ncs.array_range())
          for (const auto &cmp : crds.array_range())
            edgeCrd.push_back(cmp.as<double>());
        jsoncons::json nce = bc["Params"]["End"];
        for (const auto &crds : nce.array_range())
          for (const auto &cmp : crds.array_range())
            edgeCrd.push_back(cmp.as<double>());
        ms->FindPntsOnEdge(edgeCrd, pntIds);
        std::cout << "Number of points residing on the boundary " << bcName
                  << " is " << pntIds.size() << std::endl;
      } else {
        std::cerr << "Warning: unsupported boundary type " << bcTyp
                  << std::endl;
      }

      // register node set in Exodus II database
      if (!pntIds.empty()) {
        std::vector<int> nv;
        std::copy(pntIds.begin(), pntIds.end(), std::back_inserter(nv));
        em->addNdeSetByNdeIdLst(bcName, nv);
      }
    }
  } else if (opr == "Merge Nodes") {
    em->mergeNodes();
  } else if (opr == "Mesh Scaling") {
    em->scaleNodes(ppJson["Scale"].as<double>());
  } else {
    std::cerr << "Unknown operation requested: " << opr << std::endl;
  }
}

jsoncons::string_view ConversionDriver::getProgramType() const {
  return programType;
}

}  // namespace DRV
}  // namespace NEM
