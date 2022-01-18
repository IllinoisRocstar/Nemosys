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
#ifndef NEMOSYS_VTKMESH_H_
#define NEMOSYS_VTKMESH_H_

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

#include <vtkUnstructuredGrid.h>

#ifndef SWIG
// read files from any vtk format (except legacy .vtk)
template <class TReader>
vtkDataSet *ReadAnXMLOrSTLFile(const std::string &fileName) {
  vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName.c_str());
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

// read an ASCII legacy vtk file containing an unstructured grid
// the available vtk readers do not suffice for this purpose if the file
// contains more than one point or cell data array. This function has mainly
// been tested on vtk output from MFEM
vtkSmartPointer<vtkUnstructuredGrid> ReadALegacyVTKFile(
    const std::string &fileName);

vtkSmartPointer<vtkUnstructuredGrid> ReadDegenerateVTKFile(
    const std::string &fileName);

// helpers for reading legacy vtk

// helper that reads/checks header from legacy vtk
bool readLegacyVTKHeader(const std::string &line);
// helper that reads field data from legacy vtk
bool readLegacyVTKFieldData(const std::istream &meshStream,
                            const std::string &line,
                            vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that reads points from legacy vtk
bool readLegacyVTKPoints(const std::istream &meshStream,
                         const std::string &line, nemId_t &numPoints,
                         vtkSmartPointer<vtkPoints> points,
                         vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that reads cells from legacy vtk
bool readLegacyVTKCells(
    const std::istream &meshStream, const std::string &line, nemId_t &numCells,
    const std::vector<vtkSmartPointer<vtkIdList>> &vtkCellIds,
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that reads point and cell data from legacy vtk
bool readLegacyVTKData(const std::ifstream &meshStream, const std::string &line,
                       nemId_t numTuple, bool pointOrCell, bool &hasPointOrCell,
                       vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that casts data arrays to type specified in legacy file (not used rn)
void addLegacyVTKData(vtkDataArray *arr, const std::string &type,
                      bool pointOrCell,
                      vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);

// method to write vtk grids
template <class TWriter>
void writeVTFile(const std::string &fname,
                 vtkSmartPointer<vtkDataSet> dataSet) {
  vtkSmartPointer<TWriter> Writer = vtkSmartPointer<TWriter>::New();
  Writer->SetFileName(fname.c_str());
  Writer->SetInputData(dataSet);
  // WARNING: ASCII does not work for STL and legacy VTK files
  // WARNING: large files will be generated
  // Writer->SetDataModeToAscii();
  Writer->Update();
  Writer->Write();
}
#endif  // SWIG

class NEMOSYS_DEPRECATED_EXPORT vtkMesh
    : public meshBase {
  // constructor and destructor
 public:
  vtkMesh() = default;
  explicit vtkMesh(const std::string &fname);
  // limited support for conversion between formats
  vtkMesh(const std::string &fname1, const std::string &fname2);
  // put existing vtkDataSet into meshBase
  vtkMesh(vtkSmartPointer<vtkDataSet> dataSet_tmp, const std::string &fname);
  /* create from coordinates and connectivities.
     NOTE: use of this is only valid when mesh has one cell type.
           cellType one of the vtkCellType enums.
           Currently, only VTK_TETRA and VTK_TRIANGLE are supported */
  vtkMesh(const std::vector<double> &xCrds, const std::vector<double> &yCrds,
          const std::vector<double> &zCrds,
          const std::vector<nemId_t> &elemConn, int cellType,
          const std::string &newname);

  ~vtkMesh() override { std::cout << "vtkMesh destroyed" << std::endl; }

  // access
 public:
  // abstract read method reserved for derived classes
  void read(const std::string &fname) override {}
  // merge vtk mesh with current mesh
  void merge(vtkSmartPointer<vtkDataSet> dataSet_new);
  // get point with id
  std::vector<double> getPoint(nemId_t id) const override;
  // get 3 vecs with x,y and z coords
  std::vector<std::vector<double>> getVertCrds() const override;
  // get cell with id : returns point indices and respective coordinates
  std::map<nemId_t, std::vector<double>> getCell(nemId_t id) const override;
  std::vector<std::vector<double>> getCellVec(nemId_t id) const override;
  // get diameter of circumsphere of each cell
  std::vector<double> getCellLengths() const override;
  // get center of a cell
  std::vector<double> getCellCenter(nemId_t cellID) const override;
  // get cell type as an integer
  // assumes all elements are the same type
  int getCellType() const override;
  // get edge lengths of dataSet
  void inspectEdges(const std::string &ofname) const override;
  std::vector<nemId_t> getConnectivities() const override;

  // processing
 public:
  vtkSmartPointer<vtkDataSet> extractSurface() override;

  // diagnostics
 public:
  void report() const override;
  void write() const override { meshBase::write(); }
  void write(const std::string &fname) const override;

  // set and get point and cell data
 public:
  // set point data (numComponents per point determined by dim of data[0]
  void setPointDataArray(const std::string &name,
                         const std::vector<std::vector<double>> &data) override;
  void setPointDataArray(const std::string &name,
                         const std::vector<double> &data) override;
  // set cell data (numComponents per cell determined by dim of data[0])
  void setCellDataArray(const std::string &name,
                        const std::vector<std::vector<double>> &data) override;
  // set scalar cell data
  void setCellDataArray(const std::string &name,
                        const std::vector<double> &data) override;
  // remove point data with given id from dataSet if it exists
  void unsetPointDataArray(int arrayID) override;
  void unsetPointDataArray(const std::string &name) override;
  // remove cell data with given id from dataSet if it exists
  void unsetCellDataArray(int arrayID) override;
  void unsetCellDataArray(const std::string &name) override;
  // remove field data with given id from dataSet
  void unsetFieldDataArray(const std::string &name) override;
  // get scalar point or cell data array. assumes data is not allocated prior to
  // calling
  void getPointDataArray(const std::string &name, std::vector<double> &data)
      override;
  void getPointDataArray(int arrayId, std::vector<double> &data) override;
  int getCellDataIdx(const std::string &name) override;
  void getCellDataArray(const std::string &name, std::vector<double> &data)
      override;
  void getCellDataArray(int arrayId, std::vector<double> &data) override;
};

#endif  // NEMOSYS_VTKMESH_H_
