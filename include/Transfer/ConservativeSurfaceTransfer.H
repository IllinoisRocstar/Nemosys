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
#ifndef NEMOSYS_CONSERVATIVESURFACETRANSFER_H_
#define NEMOSYS_CONSERVATIVESURFACETRANSFER_H_

#include "Transfer/TransferBase.H"
#include "Mesh/meshBase.H"

class NEMOSYS_EXPORT ConservativeSurfaceTransfer : public TransferBase {
public:
  ConservativeSurfaceTransfer(meshBase *_source, meshBase *_target);

  ~ConservativeSurfaceTransfer() override {
    std::cout << "Conservative Surface Transfer destroyed" << std::endl;
  }

public:
  static ConservativeSurfaceTransfer *Create(meshBase *_source,
                                             meshBase *_target) {
    return new ConservativeSurfaceTransfer(_source, _target);
  }

  static std::shared_ptr<ConservativeSurfaceTransfer>
  CreateShared(meshBase *_source, meshBase *_target) {
    return std::shared_ptr<ConservativeSurfaceTransfer>(
        ConservativeSurfaceTransfer::Create(_source, _target));
  }

public:
  int transferPointData(
      const std::vector<int> &arrayIDs = std::vector<int>(),
      const std::vector<std::string> &newnames = std::vector<std::string>());

  int transferPointData(int i, vtkSmartPointer<vtkGenericCell> genCell,
                        std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSource,
                        std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget,
                        bool flip) {
    return 0;
  }

  int writeOverlay();

public:
  int transferCellData(
      const std::vector<int> &arrayIDs,
      const std::vector<std::string> &newnames = std::vector<std::string>()) {
    return 0;
  }

  int transferCellData(
      int i, vtkSmartPointer<vtkGenericCell> genCell,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSourceToPoint,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget) {
    return 0;
  }

  int run(
      const std::vector<std::string> &newnames = std::vector<std::string>()) {
    return 0;
  }

  void transfer(int arrayId);

private:

  void extractDataFromVTK(vtkDataSet *data, std::vector<double> &coords,
                          std::vector<int> &elems);

  // function handles for COM (IMPACT)
  int RFC_clear, RFC_read, RFC_write, RFC_overlay, RFC_transfer, RFC_interp;

  /*
  vtkDataSet *sourceDataSet;
  vtkDataSet *targetDataSet;

  // mesh data in COM format
  std::vector<double> sourceCoords;
  std::vector<int> sourceElems;

  std::vector<double> targetCoords;
  std::vector<int> targetElems;
  */

};

#endif  // NEMOSYS_CONSERVATIVESURFACETRANSFER_H_
