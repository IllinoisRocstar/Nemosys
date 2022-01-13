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
#ifndef NEMOSYS_MERGECELLS_H_
#define NEMOSYS_MERGECELLS_H_

#include "nemosys_export.h"

#include <map>

#include <vtkObject.h>
#include <vtkDataSetAttributes.h>
#include <vtkUnstructuredGrid.h>

class vtkDataSet;
class vtkPointData;
class vtkCellData;

namespace NEM {
namespace MSH {

class mergeCellsSTLCloak;

/*
 * Adopted from vtkMergeCells, VTK, with very minor
 * modifications to return a map between old and new cells. Note duplicating
 * the entire class becomes unnecessary as of VTK version 9, where most private
 * members of this class are protected.
 */
class NEMOSYS_NO_EXPORT mergeCells : public vtkObject {
 public:
  vtkTypeMacro(mergeCells, vtkObject)
  void PrintSelf(ostream &os, vtkIndent indent) override;

  static mergeCells *New();

  vtkSetObjectMacro(UnstructuredGrid, vtkUnstructuredGrid);
  vtkGetObjectMacro(UnstructuredGrid, vtkUnstructuredGrid);

  vtkSetMacro(TotalNumberOfCells, vtkIdType);
  vtkGetMacro(TotalNumberOfCells, vtkIdType);

  vtkSetMacro(TotalNumberOfPoints, vtkIdType);
  vtkGetMacro(TotalNumberOfPoints, vtkIdType);

  vtkSetMacro(UseGlobalIds, int);
  vtkGetMacro(UseGlobalIds, int);

  vtkSetClampMacro(PointMergeTolerance, float, 0.0f, VTK_FLOAT_MAX);
  vtkGetMacro(PointMergeTolerance, float);

  vtkSetMacro(UseGlobalCellIds, int);
  vtkGetMacro(UseGlobalCellIds, int);

  vtkSetMacro(MergeDuplicatePoints, int);
  vtkGetMacro(MergeDuplicatePoints, int);
  vtkBooleanMacro(MergeDuplicatePoints, int);

  vtkSetMacro(TotalNumberOfDataSets, int);
  vtkGetMacro(TotalNumberOfDataSets, int);

  /**
   * Provide a DataSet to be merged in to the final UnstructuredGrid.
   * This call returns after the merge has completed.  Be sure to call
   * SetTotalNumberOfCells, SetTotalNumberOfPoints, and SetTotalNumberOfDataSets
   * before making this call. Return the idMap from old and new points.
   */

  vtkSmartPointer<vtkIdTypeArray> MergeDataSet(vtkDataSet *set);

  void Finish();

 protected:
  mergeCells();
  ~mergeCells() VTK_OVERRIDE;

  void FreeLists();
  void StartUGrid(vtkDataSet *set);
  vtkIdTypeArray *MapPointsToIdsUsingGlobalIds(vtkDataSet *set);
  vtkIdTypeArray *MapPointsToIdsUsingLocator(vtkDataSet *set);
  vtkIdType AddNewCellsUnstructuredGrid(vtkDataSet *set, vtkIdTypeArray *idMap);
  vtkIdType AddNewCellsDataSet(vtkDataSet *set, vtkIdTypeArray *idMap);

  vtkIdType GlobalCellIdAccessGetId(vtkIdType idx);
  int GlobalCellIdAccessStart(vtkDataSet *set);
  vtkIdType GlobalNodeIdAccessGetId(vtkIdType idx);
  int GlobalNodeIdAccessStart(vtkDataSet *set);

  int TotalNumberOfDataSets;

  vtkIdType TotalNumberOfCells;
  vtkIdType TotalNumberOfPoints;

  vtkIdType NumberOfCells;  // so far
  vtkIdType NumberOfPoints;

  int UseGlobalIds;  // point, or node, IDs
  int GlobalIdArrayType;
  void *GlobalIdArray;

  int UseGlobalCellIds;  // cell IDs
  int GlobalCellIdArrayType;
  void *GlobalCellIdArray;

  float PointMergeTolerance;
  int MergeDuplicatePoints;

  char InputIsUGrid;
  char InputIsPointSet;

  mergeCellsSTLCloak *GlobalIdMap;
  mergeCellsSTLCloak *GlobalCellIdMap;

  vtkDataSetAttributes::FieldList *ptList;
  vtkDataSetAttributes::FieldList *cellList;

  vtkUnstructuredGrid *UnstructuredGrid;

  int nextGrid;

  mergeCells(const mergeCells &) = delete;
  void operator=(const mergeCells &) = delete;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_MERGECELLS_H_
