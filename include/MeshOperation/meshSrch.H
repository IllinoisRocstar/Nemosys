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
#ifndef NEMOSYS_MESHSRCH_H_
#define NEMOSYS_MESHSRCH_H_

#include <set>

#include <vtkCellLocator.h>

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

class NEMOSYS_EXPORT meshSrch : public meshBase {
  // constructors and destructors
 public:
  meshSrch() = delete;
  explicit meshSrch(meshBase *mb) : meshBase(*mb), upd_vcl(true) {}

  static meshSrch *Create(meshBase *mb) {
    auto *ms = new meshSrch(mb);
    return ms;
  }

  static std::shared_ptr<meshSrch> CreateShared(meshBase *mb) {
    std::shared_ptr<meshSrch> ms;
    ms.reset(meshSrch::Create(mb));
    return ms;
  }

  static std::unique_ptr<meshSrch> CreateUnique(meshBase *mb) {
    return std::unique_ptr<meshSrch>(meshSrch::Create(mb));
  }

  ~meshSrch() override = default;

  // general access
 public:
  std::vector<std::vector<double>> getCellVec(nemId_t id) const override;
  std::vector<double> getPoint(nemId_t id) const override;
  std::vector<double> getCellCenter(nemId_t cellID) const override;

  int getCellType() const override { return 0; }
  std::vector<std::vector<double>> getVertCrds() const override {
    return std::vector<std::vector<double>>();
  }
  std::map<nemId_t, std::vector<double>> getCell(nemId_t id) const override {
    return std::map<nemId_t, std::vector<double>>();
  }
  int getCellDataIdx(const std::string &name) override { return 0; }
  std::vector<double> getCellLengths() const override {
    return std::vector<double>();
  }
  std::vector<nemId_t> getConnectivities() const override {
    return std::vector<nemId_t>();
  }

  // check for special conditions
 public:
  bool chkDuplElm() const;  // finds duplicate elements

  // point search methods
 public:
  // get coordinates and connectivities of the surface triangulation and returns
  // ids for the nodes that reside on the triangulation within given tolerance.
  void FindPntsOnTriSrf(const std::vector<double> &crds,
                        const std::vector<nemId_t> &conn,
                        std::set<nemId_t> &ids, double tol = 0.1e-15) const;
  // get coordinates of the start and end points of an edge and returns ids for
  // the nodes that reside on the edge within given tolerance.
  void FindPntsOnEdge(std::vector<double> &crds, std::set<nemId_t> &ids,
                      double tol = 0.1e-15) const;

  // cell search methods
 public:
  // get vtkPolyData and return ids for the cells whose center reside inside the
  // polyData.
  void FindCellsInPolyData(vtkPolyData *polyData, std::vector<nemId_t> &ids,
                           bool query3Donly = true, double tol = 0.1e-15) const;
  // get coordinates of a bounding box [xmin,xmax, ymin,ymax, zmin,zmax] and
  // return ids for the cells whose center reside inside the bounding box.
  void FindCellsWithinBounds(std::vector<double> &bb, std::vector<nemId_t> &ids,
                             bool fulImrsd = true);
  // get coordinates and connectivities of the surface triangulation of a closed
  // manifold and return ids for the cells whose center reside inside the
  // manifold.
  void FindCellsInTriSrf(const std::vector<std::vector<double>> &crds,
                         const std::vector<std::vector<vtkIdType>> &conns,
                         std::vector<nemId_t> &ids, bool query3Donly = true,
                         double tol = 0.1e-15) const;
  // get center and radius of sphere and return ids for the cells whose center
  // reside inside the sphere.
  void FindCellsInSphere(const std::vector<double> &center, double radius,
                         std::vector<nemId_t> &ids, bool query3Donly = true,
                         double tol = 0.1e-15) const;

  // misc
 public:
  void inspectEdges(const std::string &ofname) const override {}
  vtkSmartPointer<vtkDataSet> extractSurface() override {
    return vtkSmartPointer<vtkDataSet>();
  }

  void read(const std::string &fname = std::string()) override {}
  using meshBase::write;
  void write(const std::string &fname) const override {}

  // internal management
 private:
  void buildCellLocator();

 private:
  bool upd_vcl;
  vtkSmartPointer<vtkCellLocator> vcl;
};

#endif  // NEMOSYS_MESHSRCH_H_
