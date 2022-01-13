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
#ifndef NEMOSYS_CUBATURE_H_
#define NEMOSYS_CUBATURE_H_

#include <memory>

// Nemosys headers
#include "nemosys_export.h"

// VTK
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkQuadratureSchemeDefinition.h>

// 3 point symmetric gaussian quadrature for triangle element
extern double TRI3[];
extern double TRI3W[];
// 4 point symmetric gaussian quadrature for tetrahedral element
extern double TET4[];
extern double TET4W[];

// pair type for coordinate and data there
using pntDataPair = std::pair<std::vector<double>, std::vector<double>>;
// holds gauss points and data at these points as pairs
using pntDataPairVec = std::vector<pntDataPair>;

class NEMOSYS_EXPORT GaussCubature {
 public:
  explicit GaussCubature(vtkDataSet *_dataSet);
  GaussCubature(vtkDataSet *_dataSet, const std::vector<int> &arrayIDs);

  ~GaussCubature() { delete[] dict; }

  GaussCubature(const GaussCubature &that) = delete;
  GaussCubature &operator=(const GaussCubature &that) = delete;

  void constructGaussMesh();

 public:
  /**
   * @brief Returns coordinates of gauss points and associated data at cell
   * @param cellId cell id
   * @return pntDataPairVec (see alias in Cubature.H)
   */
  pntDataPairVec getGaussPointsAndDataAtCell(int cellID);
  /**
   * @brief Interpolate values to gauss points for arrays specified in arrayIDs
   * at construction time.
   */
  void interpolateToGaussPoints();
  /**
   * @brief Interpolate values to gauss points for arrays specified in arrayIDs
   * at construction time and specify new names (in order).
   * @param newArrayNames name of transferred data from arrayIDs[i] will be
   * named newArrayNames[i]
   */
  void interpolateToGaussPoints(const std::vector<std::string> &newArrayNames);
  /**
   * @brief Integrate arrays specified by arrayIDs over all cells.
   * @return Integral data (first index - array, second index - component)
   */
  std::vector<std::vector<double>> integrateOverAllCells();
  /**
   * @brief Same as integrateOverAllCells() but with integral value array
   * names specified by newArrayNames and an optional flag to compute
   * the root mean square error (RMSE) integral over an array with error values,
   * an error metric normalized by cell volume.
   * @param newArrayNames integral value array names
   * @param computeRMSE compute RMSE (see above)
   * @return Integral data (first index - array, second index - component)
   */
  std::vector<std::vector<double>> integrateOverAllCells(
      const std::vector<std::string> &newArrayNames, bool computeRMSE);

  void setArrayIDs(const std::vector<int> &_arrayIDs) { arrayIDs = _arrayIDs; }
  double computeJacobian(vtkSmartPointer<vtkGenericCell> genCell,
                         int cellType) const;
  double computeCellVolume(vtkSmartPointer<vtkGenericCell> genCell,
                           int cellType) const;

  // access
 public:
  vtkDataSet *getDataSet() const { return dataSet; }

  vtkSmartPointer<vtkPolyData> getGaussMesh() const { return gaussMesh; }
  vtkQuadratureSchemeDefinition **getDict() const { return dict; }
  std::vector<int> getNumComponents() const { return numComponents; }
  std::vector<int> getArrayIDs() const { return arrayIDs; }
  int getTotalComponents() const { return totalComponents; }

  void writeGaussMesh(const char *name) const;

  // factory constructors
 public:
  static GaussCubature *Create(vtkDataSet *_dataSet);
  static GaussCubature *Create(vtkDataSet *_dataSet,
                               const std::vector<int> &arrayIDs);
  static std::unique_ptr<GaussCubature> CreateUnique(vtkDataSet *_dataSet);
  static std::unique_ptr<GaussCubature> CreateUnique(
      vtkDataSet *_dataSet, const std::vector<int> &arrayIDs);
  static std::shared_ptr<GaussCubature> CreateShared(vtkDataSet *_dataSet);
  static std::shared_ptr<GaussCubature> CreateShared(
      vtkDataSet *_dataSet, const std::vector<int> &arrayIDs);

 private:
  // we want gauss points of this mesh
  vtkSmartPointer<vtkDataSet> dataSet;
  // number of volume elements in node mesh
  int numVolCells;
  // we put the gauss points and interpolated data into this mesh
  vtkSmartPointer<vtkPolyData> gaussMesh;
  // dictionary relating cell type quadrature scheme info for that type
  vtkQuadratureSchemeDefinition **dict;
  // array ids of data to be interpolated
  std::vector<int> arrayIDs;
  std::vector<int> numComponents;
  int totalComponents;

  // get offset from nodeMesh for lookup of gauss points in gaussMesh
  int getOffset(int cellID) const;
  // interpolates provided data (das) to gauss points in cell
  int interpolateToGaussPointsAtCell(
      int cellID, vtkSmartPointer<vtkGenericCell> genCell,
      const std::vector<vtkSmartPointer<vtkDataArray>> &das,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &daGausses) const;
  // integrates provided data over cell
  void integrateOverCell(
      int cellID, vtkSmartPointer<vtkGenericCell> genCell,
      vtkSmartPointer<vtkPointData> pd,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &integralData,
      std::vector<std::vector<double>> &totalIntegralData) const;
  // integrates data specified by newArrayNames over cell,
  // normalizes the integral by the cell volume and takes the sqrt
  void integrateOverCell(
      int cellID, vtkSmartPointer<vtkGenericCell> genCell,
      vtkSmartPointer<vtkPointData> pd,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &integralData,
      std::vector<std::vector<double>> &totalIntegralData,
      const std::vector<std::string> &newArrayNames, bool normalizeByVol) const;
};

#endif  // NEMOSYS_CUBATURE_H_
