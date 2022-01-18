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
#include "MeshOperation/meshSrch.H"

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellCenters.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSphereSource.h>

// DEBUG:
//#include <vtkSTLWriter.h>

#include "AuxiliaryFunctions.H"

using nemAux::operator*;  // for vector multiplication.
using nemAux::operator+;  // for vector addition.

// get point with id
std::vector<double> meshSrch::getPoint(nemId_t id) const {
  double coords[3];
  dataSet->GetPoint(id, coords);
  std::vector<double> result(coords, coords + 3);
  return result;
}

// returns coordinates of the cell vertices in a vector
std::vector<std::vector<double>> meshSrch::getCellVec(nemId_t id) const {
  if (id < numCells) {
    std::vector<std::vector<double>> cell;
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    point_ids = dataSet->GetCell(id)->GetPointIds();
    vtkIdType num_ids = point_ids->GetNumberOfIds();
    cell.resize(num_ids);
    for (vtkIdType i = 0; i < num_ids; ++i) {
      nemId_t pntId = point_ids->GetId(i);
      cell[i] = getPoint(pntId);
    }
    return cell;
  } else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(1);
  }
}

// get center of a cell
std::vector<double> meshSrch::getCellCenter(nemId_t cellID) const {
  std::vector<double> center(3);
  std::vector<std::vector<double>> cell = getCellVec(cellID);

  for (const auto &i : cell) center = center + i;
  return (1.0 / static_cast<double>(cell.size())) * center;
}

void meshSrch::buildCellLocator() {
  if (upd_vcl) {
    // Create the tree
    vcl = vtkSmartPointer<vtkCellLocator>::New();
    vcl->SetDataSet(dataSet);
    vcl->BuildLocator();
    upd_vcl = false;
  }
}

void meshSrch::FindCellsWithinBounds(std::vector<double> &bb,
                                     std::vector<nemId_t> &ids, bool fulImrsd) {
  // finding all intersecting cells
  buildCellLocator();
  vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
  vcl->FindCellsWithinBounds(bb.data(), idl);
  std::cout << "Found " << idl->GetNumberOfIds() << " cells." << std::endl;
  std::vector<nemId_t> aids;
  for (vtkIdType idx = 0; idx < idl->GetNumberOfIds(); idx++)
    aids.push_back(idl->GetId(idx));
  // removing cells centered out of the bounding box
  int nr = 0;
  if (fulImrsd) {
    for (const auto &aid : aids)
      if (!nemAux::isInBBox(getCellCenter(aid), bb)) {
        nr++;
        continue;
      } else {
        ids.push_back(aid);
      }
    std::cout << "Remove " << nr << " cells from the list." << std::endl;
  } else {
    ids = aids;
  }
}

void meshSrch::FindPntsOnTriSrf(const std::vector<double> &crds,
                                const std::vector<nemId_t> &conn,
                                std::set<nemId_t> &ids, double tol) const {
  // create polyData
  vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
  for (std::size_t iPnt = 0; iPnt < crds.size() / 3; iPnt++)
    pnts->InsertNextPoint(crds[iPnt * 3], crds[iPnt * 3 + 1],
                          crds[iPnt * 3 + 2]);
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  for (std::size_t iCel = 0; iCel < conn.size() / 3; iCel++) {
    polys->InsertNextCell(3);
    polys->InsertCellPoint(conn[iCel * 3]);
    polys->InsertCellPoint(conn[iCel * 3 + 1]);
    polys->InsertCellPoint(conn[iCel * 3 + 2]);
  }
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(pnts);
  polyData->SetPolys(polys);

  // Write the file
  // vtkSmartPointer<vtkXMLPolyDataWriter> writer
  //    = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  // writer->SetFileName("test.vtp");
  // writer->SetInputData(polyData);
  // Optional - set the mode. The default is binary.
  ////writer->SetDataModeToBinary();
  // writer->SetDataModeToAscii();
  // writer->Write();

  // find nodes residing on the trisurf
  // create cell locator
  vtkSmartPointer<vtkCellLocator> cellLocator =
      vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(polyData);
  cellLocator->BuildLocator();

  double closestPoint[3];
  double closestPointDist2;
  vtkIdType cellId;
  int subId;
  for (nemId_t iPt = 0; iPt < getNumberOfPoints(); iPt++) {
    std::vector<double> pnt = getPoint(iPt);
    cellLocator->FindClosestPoint(pnt.data(), closestPoint, cellId, subId,
                                  closestPointDist2);
    if (closestPointDist2 < tol) ids.insert(iPt + 1);
  }
}

void meshSrch::FindPntsOnEdge(std::vector<double> &crds, std::set<nemId_t> &ids,
                              double tol) const {
  // create polyData
  vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
  for (std::size_t iPnt = 0; iPnt < crds.size() / 3; iPnt++)
    pnts->InsertNextPoint(crds[iPnt * 3], crds[iPnt * 3 + 1],
                          crds[iPnt * 3 + 2]);
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  polys->InsertNextCell(2);
  polys->InsertCellPoint(0);
  polys->InsertCellPoint(1);
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(pnts);
  polyData->SetPolys(polys);

  // Write the file
  // vtkSmartPointer<vtkXMLPolyDataWriter> writer
  //    = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  // writer->SetFileName("edge.vtp");
  // writer->SetInputData(polyData);
  // Optional - set the mode. The default is binary.
  ////writer->SetDataModeToBinary();
  // writer->SetDataModeToAscii();
  // writer->Write();

  // find nodes residing on the edge
  // create cell locator
  vtkSmartPointer<vtkCellLocator> cellLocator =
      vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(polyData);
  cellLocator->BuildLocator();

  double closestPoint[3];
  double closestPointDist2;
  vtkIdType cellId;
  int subId;
  for (nemId_t iPt = 0; iPt < getNumberOfPoints(); iPt++) {
    std::vector<double> pnt = getPoint(iPt);
    cellLocator->FindClosestPoint(pnt.data(), closestPoint, cellId, subId,
                                  closestPointDist2);
    if (closestPointDist2 < tol) ids.insert(iPt + 1);
  }
}

// checks for duplicate elements
bool meshSrch::chkDuplElm() const {
  std::set<std::vector<nemId_t>> ids;
  for (nemId_t ic = 0; ic < getNumberOfCells(); ic++) {
    std::vector<nemId_t> cid;
    vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
    idl = dataSet->GetCell(ic)->GetPointIds();
    for (vtkIdType id = 0; id < idl->GetNumberOfIds(); id++)
      cid.push_back(idl->GetId(id));
    std::pair<std::set<std::vector<nemId_t>>::iterator, bool> ret;
    ret = ids.insert(cid);
    if (!ret.second) return true;
  }
  return false;
}

void meshSrch::FindCellsInPolyData(vtkPolyData *polyData,
                                   std::vector<nemId_t> &ids, bool query3Donly,
                                   double tol) const {
  // pass dataSet through vtkCellCenters filter
  vtkSmartPointer<vtkCellCenters> cc = vtkSmartPointer<vtkCellCenters>::New();

  cc->SetInputData(dataSet);

  // create vtkSelectEnclosedPoints
  vtkSmartPointer<vtkSelectEnclosedPoints> sep =
      vtkSmartPointer<vtkSelectEnclosedPoints>::New();

  sep->SetInputConnection(cc->GetOutputPort());

  sep->SetSurfaceData(polyData);
  // sep->CheckSurfaceOn();

  sep->SetTolerance(tol);

  sep->Update();

  for (nemId_t id = 0; id < getNumberOfCells(); ++id)
    if (sep->IsInside(id)) {
      if (query3Donly && dataSet->GetCell(id)->GetCellDimension() != 3)
        continue;
      ids.emplace_back(id);
    }
}

void meshSrch::FindCellsInTriSrf(
    const std::vector<std::vector<double>> &crds,
    const std::vector<std::vector<vtkIdType>> &conns, std::vector<nemId_t> &ids,
    bool query3Donly, double tol) const {
  // create vtkPolyData using crds and conns
  vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

  for (const auto &crd : crds) pnts->InsertNextPoint(crd.data());
  for (const auto &conn : conns)
    polys->InsertNextCell(conn.size(), conn.data());
  polyData->SetPoints(pnts);
  polyData->SetPolys(polys);

  // DEBUG: Write polyData to STL
  //  vtkSmartPointer<vtkSTLWriter> stlw = vtkSmartPointer<vtkSTLWriter>::New();
  //  stlw->SetInputData(polyData);
  //  stlw->SetFileName("polydata.stl");
  //  stlw->Write();

  FindCellsInPolyData(polyData, ids, query3Donly, tol);
}

void meshSrch::FindCellsInSphere(const std::vector<double> &center,
                                 double radius, std::vector<nemId_t> &ids,
                                 bool query3Donly, double tol) const {
  // create vtkSphere using center and radius
  vtkSmartPointer<vtkSphereSource> ss = vtkSmartPointer<vtkSphereSource>::New();
  ss->SetRadius(radius);
  ss->SetCenter(center[0], center[1], center[2]);

  ss->Update();

  FindCellsInPolyData(ss->GetOutput(), ids, query3Donly, tol);
}
