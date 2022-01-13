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
#include "Transfer/ConservativeSurfaceTransfer.H"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <com.h>

#include "AuxiliaryFunctions.H"

// debug
#include <iostream>

COM_EXTERN_MODULE(SurfX)
COM_EXTERN_MODULE(SimOut)

ConservativeSurfaceTransfer::ConservativeSurfaceTransfer(meshBase *_source,
                                                         meshBase *_target) {
  source = _source;
  target = _target;
}

int ConservativeSurfaceTransfer::transferPointData(
    const std::vector<int> &arrayIds,
    const std::vector<std::string> &newnames) {

  COM_init(NULL, NULL);
  COM_set_profiling(1);
  // COM_set_verbose(10);

  MPI_Comm comm = MPI_COMM_WORLD;

  COM_LOAD_MODULE_STATIC_DYNAMIC(SurfX, "RFC");

  RFC_clear = COM_get_function_handle("RFC.clear_overlay");
  RFC_read = COM_get_function_handle("RFC.read_overlay");
  RFC_write = COM_get_function_handle("RFC.write_overlay");
  RFC_overlay = COM_get_function_handle("RFC.overlay");
  RFC_transfer = COM_get_function_handle("RFC.least_squares_transfer");
  RFC_interp = COM_get_function_handle("RFC.interpolate");

  vtkDataSet *sourceDataSet = source->getDataSet();
  vtkDataSet *targetDataSet = target->getDataSet();

  std::vector<double> sourceCoords;
  std::vector<int> sourceElems;
  extractDataFromVTK(sourceDataSet, sourceCoords, sourceElems);

  // === SOURCE WINDOW ===

  COM_new_window("sourceWindow");

  // set source nodal coordinates (nc)
  COM_set_size("sourceWindow.nc", 1, sourceCoords.size() / 3);
  COM_set_array("sourceWindow.nc", 1, &sourceCoords[0], 3);
  // set source element node ids (:t3:)
  COM_set_size("sourceWindow.:t3:", 1, sourceElems.size() / 3);
  COM_set_array("sourceWindow.:t3:", 1, &sourceElems[0], 3);

  // initialize source soln array - this is reused
  COM_new_dataitem("sourceWindow.soln", 'n', COM_DOUBLE, 1, "");

  COM_window_init_done("sourceWindow");

  // === TARGET WINDOW ===

  std::vector<double> targetCoords;
  std::vector<int> targetElems;
  extractDataFromVTK(target->getDataSet(), targetCoords, targetElems);

  COM_new_window("targetWindow");
  // set target nodal coordinates (nc)
  COM_set_size("targetWindow.nc", 1, targetCoords.size() / 3);
  COM_set_array("targetWindow.nc", 1, &targetCoords[0], 3);
  // set target element node ids (:t3:)
  COM_set_size("targetWindow.:t3:", 1, targetElems.size() / 3);
  COM_set_array("targetWindow.:t3:", 1, &targetElems[0], 3);

  // initialize target soln array - this is reused
  std::vector<double> targetData(targetCoords.size(), -1.);
  COM_new_dataitem("targetWindow.soln", 'n', COM_DOUBLE, 1, "m/s");
  COM_resize_array("targetWindow.soln");

  COM_window_init_done("targetWindow");

  // === OVERLAY ===
  int srcMesh = COM_get_dataitem_handle("sourceWindow.mesh");
  int tgtMesh = COM_get_dataitem_handle("targetWindow.mesh");
  COM_call_function(RFC_overlay, &srcMesh, &tgtMesh);

  // === TRANSFER ===
  // set source data and transfer
  int srcData = COM_get_dataitem_handle("sourceWindow.soln");
  int tgtData = COM_get_dataitem_handle("targetWindow.soln");

  // source data buffer
  std::vector<double> sourceData(sourceCoords.size());
  // transfer specified arrays
  for (const int &arrayId : arrayIds) {
    vtkSmartPointer<vtkDataArray> sourceArray =
        sourceDataSet->GetPointData()->GetArray(arrayId);

    auto targetPointValues = vtkSmartPointer<vtkDoubleArray>::New();
    targetPointValues->SetName(sourceArray->GetName());
    targetPointValues->SetNumberOfComponents(
        sourceArray->GetNumberOfComponents());

    for (int componentId = 0;
         componentId < sourceArray->GetNumberOfComponents(); ++componentId) {
      auto sourcePointValues = vtkDoubleArray::SafeDownCast(
          sourceDataSet->GetPointData()->GetArray(arrayId));

      for (int srcNodeId = 0; srcNodeId < sourceCoords.size() / 3;
           ++srcNodeId) {
        sourceData[srcNodeId] =
            sourcePointValues->GetComponent(srcNodeId, componentId);
      }

      COM_set_array("sourceWindow.soln", 1, &sourceData[0]);

      // transfer
      COM_call_function(RFC_transfer, &srcData, &tgtData);

      void *addr;
      double *data;

      // window var, pane number, addr
      COM_get_array("targetWindow.soln", 1, &addr);
      data = (double *)addr;

      for (int trgNodeId = 0; trgNodeId < targetCoords.size() / 3;
           ++trgNodeId) {
        targetPointValues->InsertComponent(trgNodeId, componentId,
                                           data[trgNodeId]);
      }
    }

    targetDataSet->GetPointData()->AddArray(targetPointValues);
  }

  COM_finalize();

  return 0;
}

int ConservativeSurfaceTransfer::writeOverlay()
{
  // TODO Find a way to resolve scoping with COM
  // this initialization should not need to happen in more than one function
  COM_init(NULL, NULL);
  COM_set_profiling(1);
  // COM_set_verbose(10);

  MPI_Comm comm = MPI_COMM_WORLD;

  COM_LOAD_MODULE_STATIC_DYNAMIC(SurfX, "RFC");

  RFC_write = COM_get_function_handle("RFC.write_overlay");
  RFC_overlay = COM_get_function_handle("RFC.overlay");

  vtkDataSet *sourceDataSet = source->getDataSet();
  vtkDataSet *targetDataSet = target->getDataSet();

  std::vector<double> sourceCoords;
  std::vector<int> sourceElems;
  extractDataFromVTK(sourceDataSet, sourceCoords, sourceElems);

  // === SOURCE WINDOW ===
  COM_new_window("sourceWindow");

  // set source nodal coordinates (nc)
  COM_set_size("sourceWindow.nc", 1, sourceCoords.size() / 3);
  COM_set_array("sourceWindow.nc", 1, &sourceCoords[0], 3);
  // set source element node ids (:t3:)
  COM_set_size("sourceWindow.:t3:", 1, sourceElems.size() / 3);
  COM_set_array("sourceWindow.:t3:", 1, &sourceElems[0], 3);

  // initialize source soln array - this is reused
  COM_new_dataitem("sourceWindow.soln", 'n', COM_DOUBLE, 1, "");

  COM_window_init_done("sourceWindow");

  // === TARGET WINDOW ===
  std::vector<double> targetCoords;
  std::vector<int> targetElems;
  extractDataFromVTK(target->getDataSet(), targetCoords, targetElems);

  COM_new_window("targetWindow");
  // set target nodal coordinates (nc)
  COM_set_size("targetWindow.nc", 1, targetCoords.size() / 3);
  COM_set_array("targetWindow.nc", 1, &targetCoords[0], 3);
  // set target element node ids (:t3:)
  COM_set_size("targetWindow.:t3:", 1, targetElems.size() / 3);
  COM_set_array("targetWindow.:t3:", 1, &targetElems[0], 3);

  // initialize target soln array - this is reused
  std::vector<double> targetData(targetCoords.size(), -1.);
  COM_new_dataitem("targetWindow.soln", 'n', COM_DOUBLE, 1, "m/s");
  COM_resize_array("targetWindow.soln");

  COM_window_init_done("targetWindow");

  int srcMesh = COM_get_dataitem_handle("sourceWindow.mesh");
  int tgtMesh = COM_get_dataitem_handle("targetWindow.mesh");
  COM_call_function(RFC_overlay, &srcMesh, &tgtMesh);

  COM_call_function(RFC_write, &srcMesh, &tgtMesh, "srcMesh", "tgtMesh", "CGNS");

  COM_finalize();

  return 0;
}

void ConservativeSurfaceTransfer::extractDataFromVTK(
    vtkDataSet *data, std::vector<double> &coords, std::vector<int> &elems) {
  for (int ptId = 0; ptId < data->GetNumberOfPoints(); ++ptId) {
    double x[3];

    data->GetPoint(ptId, x);

    coords.push_back(x[0]);
    coords.push_back(x[1]);
    coords.push_back(x[2]);
  }

  for (int cellId = 0; cellId < data->GetNumberOfCells(); ++cellId) {
    if (data->GetCellType(cellId) != VTK_TRIANGLE)
      continue;

    vtkCell *cell = data->GetCell(cellId);

    elems.push_back(cell->GetPointId(0) + 1);
    elems.push_back(cell->GetPointId(1) + 1);
    elems.push_back(cell->GetPointId(2) + 1);
  }
}
