#ifndef NEMOSYS_GEOMESHBASE_H_
#define NEMOSYS_GEOMESHBASE_H_

#include "nemosys_export.h"

#include <array>
#include <cstdint>
#include <memory>
#include <ostream>
#include <string>

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkFieldData.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

namespace NEM {
namespace MSH {

typedef std::int64_t nemId_t;

/**
 * @class   GmshInterface
 * @brief   management class for Gmsh interface
 *
 * Gmsh uses global (and static members) variables that are initialized and
 * finalized by Gmsh API commands. This class manages geoMeshBase usage of Gmsh
 * to guarantee finalization only when no geoMeshBase instances are left.
 */
class NEMOSYS_EXPORT GmshInterface {
 public:
  ~GmshInterface();

  /**@{*/
  GmshInterface(const GmshInterface &) = delete;
  GmshInterface &operator=(const GmshInterface &) = delete;
  GmshInterface(GmshInterface &&) = delete;
  GmshInterface &operator=(GmshInterface &&) = delete;
  /**@}*/

 public:
  /**
   * Initialize Gmsh. Will only call gmsh::initialize() on first call to method.
   */
  static void Initialize();
  /**
   * Finalize Gmsh. Will only call gmsh::finalize() if count is zero.
   */
  static void Finalize();

 private:
  GmshInterface();

 private:
  static std::shared_ptr<GmshInterface> instance;
  static int count;
};

/**
 * @class   geoMeshBase
 * @brief   abstract class to specify geometry and mesh data
 *
 * geoMeshBase is an abstract class that specifies an interface for geometry
 * and mesh objects.
 */
class NEMOSYS_EXPORT geoMeshBase {
 public:
  /**
   * Convert 2D vtk dataset into gmsh. Returns gmsh tag of 2D surface added.
   *
   * Will add nodes and elements. Supports only one 0D, 1D, and 2D element set
   * at the time. Adding 3D elements will error.
   *
   * @param dataSet VTK unstructured grid
   * @return gmsh tag of 2D entity added
   */
  static int vtkSurf2gmsh(vtkUnstructuredGrid *dataSet);

  /**
   * Query for conversion from vtk to gmsh
   * @param vtkType
   * @return
   */
  static int getGmshTypeFromVTKType(int vtkType);

  /**
   * Compute discrete geometry from vtk mesh and name of geometry cell array
   * @param dataSet VTK unstructured grid
   * @param geoName name of geometry cell array
   * @param angleThreshold dihedral angle (in radians) to detect surfaces
   */
  static void computeDiscreteGeoFromMsh(vtkUnstructuredGrid *dataSet,
                                        const std::string &geoName,
                                        double angleThreshold);

 public:
  geoMeshBase();

  /**
   * Disable copy. Allow move.
   * @{
   */
  geoMeshBase(const geoMeshBase &) = delete;
  geoMeshBase &operator=(const geoMeshBase &) = delete;

  geoMeshBase(geoMeshBase &&) = default;
  geoMeshBase &operator=(geoMeshBase &&) = default;
  /**@}*/

  virtual ~geoMeshBase();

 public:
  /**
   * Write mesh to file
   * @param fileName file name
   */
  virtual void write(const std::string &fileName) = 0;

  /**
   * Print a report about the mesh
   * @param out stream
   */
  virtual void report(std::ostream &out) const = 0;
  friend std::ostream &operator<<(std::ostream &os, const geoMeshBase &base) {
    base.report(os);
    return os;
  }

 public:
  /**
   * Get the name of the geometric entities array
   * @return name of geometric entities array
   */
  const std::string &getGeoEntArrayName() const { return _geoMesh.link; }
  /**
   * Set the name of the geometric entities array
   * @param geoEntArrayName name of geometric entities array
   */
  void setGeoEntArrayName(const std::string &geoEntArrayName) {
    auto a = _geoMesh.mesh->GetCellData()->GetArray(_geoMesh.link.c_str());
    if (a) a->SetName(geoEntArrayName.c_str());
    _geoMesh.link = geoEntArrayName;
  }
  /**
   * Get the angle threshold used for surface classification
   * @return angle threshold (in radians)
   */
  double getAngleThreshold() const { return _angleThreshold; }
  /**
   * Set the angle threshold used for surface classification
   * @param angleThreshold (in radians)
   */
  void setAngleThreshold(double angleThreshold) {
    _angleThreshold = angleThreshold;
  }

 public:
  /**
   * Get the number of points in mesh
   * @return number of points in mesh
   */
  nemId_t getNumberOfPoints() const {
    return _geoMesh.mesh->GetNumberOfPoints();
  }
  /**
   * Get the number of cells in mesh
   * @return number of cells in mesh
   */
  nemId_t getNumberOfCells() const { return _geoMesh.mesh->GetNumberOfCells(); }
  /**
   * Get the coordinate of a point
   * @param id point id
   * @param x array of coordinates as (x,y,z)
   */
  void getPoint(nemId_t id, std::array<double, 3> *x) const {
    _geoMesh.mesh->GetPoint(id, x->data());
  }
  /**
   * Get the coordinate of a point
   * @param id point id
   * @return array of coordinates as (x,y,z)
   */
  std::array<double, 3> getPoint(nemId_t id) const {
    std::array<double, 3> x{};
    this->getPoint(id, &x);
    return x;
  }
  /**
   * Get cell.
   * @param cellId cell id
   * @param cell cell
   */
  void getCell(nemId_t cellId, vtkGenericCell *cell) const {
    _geoMesh.mesh->GetCell(cellId, cell);
  }
  /**
   * Get cell bounds
   * @param cellId cell id
   * @param bounds array of bounds as (x_min,x_max, y_min,y_max, z_min,z_max)
   */
  void getCellBounds(nemId_t cellId, std::array<double, 6> *bounds) const {
    _geoMesh.mesh->GetCellBounds(cellId, bounds->data());
  }
  /**
   * Get cell bounds
   * @param cellId cell id
   * @return array of bounds as (x_min,x_max, y_min,y_max, z_min,z_max)
   */
  std::array<double, 6> getCellBounds(nemId_t cellId) const {
    std::array<double, 6> bounds{};
    this->getCellBounds(cellId, &bounds);
    return bounds;
  }
  /**
   * Get VTK cell type
   * @param cellId cell id
   * @return VTK cell type (see vtk documentation for specifics)
   */
  VTKCellType getCellType(nemId_t cellId) const {
    return static_cast<VTKCellType>(_geoMesh.mesh->GetCellType(cellId));
  }
  /**
   * Get list of point ids defining specified cell
   * @param cellId cell id
   * @param ptIds list of point ids
   */
  void getCellPoints(nemId_t cellId, vtkIdList *ptIds) const {
    _geoMesh.mesh->GetCellPoints(cellId, ptIds);
  }
  /**
   * Get list of cell ids using specified point
   * @param ptId point id
   * @param cellIds list of cell ids
   */
  void getPointCells(nemId_t ptId, vtkIdList *cellIds) const {
    _geoMesh.mesh->GetPointCells(ptId, cellIds);
  }
  /**
   * TODO
   * @param cellId
   * @param ptIds
   * @param cellIds
   */
  void getCellNeighbors(nemId_t cellId, vtkIdList *ptIds,
                        vtkIdList *cellIds) const {
    _geoMesh.mesh->GetCellNeighbors(cellId, ptIds, cellIds);
  }

 public:
  /**
   * Get copy of point data array by name. Also returns array index.
   * @param arrayName array name
   * @param array point data array
   * @param arrayIdx array index
   */
  void getPointDataArrayCopy(const std::string &arrayName,
                             vtkAbstractArray *array, int *arrayIdx) const {
    array->DeepCopy(_geoMesh.mesh->GetPointData()->GetAbstractArray(
        arrayName.c_str(), *arrayIdx));
  }
  /**
   * Get copy of cell data array by name. Also returns array index.
   * @param arrayName array name
   * @param array cell data array
   * @param arrayIdx array index
   */
  void getCellDataArrayCopy(const std::string &arrayName,
                            vtkAbstractArray *array, int *arrayIdx) const {
    array->DeepCopy(_geoMesh.mesh->GetCellData()->GetAbstractArray(
        arrayName.c_str(), *arrayIdx));
  }
  /**
   * Get copy of field data array by name. Also returns array index.
   * @param arrayName array name
   * @param array field data array
   * @param arrayIdx array index
   */
  void getFieldDataArrayCopy(const std::string &arrayName,
                             vtkAbstractArray *array, int *arrayIdx) const {
    array->DeepCopy(_geoMesh.mesh->GetFieldData()->GetAbstractArray(
        arrayName.c_str(), *arrayIdx));
  }

 protected:
  struct GeoMesh {
    vtkSmartPointer<vtkUnstructuredGrid> mesh;  ///< mesh
    std::string geo;                            ///< name of gmsh model with geo

    std::string link;  ///< name of vtk int array storing geometric entities
  };

 protected:
  /**
   * Construct geoMeshBase from a GeoMesh struct
   * @param inGeoMesh a geo and mesh pair
   */
  explicit geoMeshBase(GeoMesh inGeoMesh);

 protected:
  // virtual void update();

  // void addPoints();

  /**
   * Set point data
   * @param array array pointer
   * @return array index
   */
  int setPointDataArray(vtkAbstractArray *array) {
    return _geoMesh.mesh->GetPointData()->AddArray(array);
  }
  /**
   * Set cell data
   * @param array array pointer
   * @return array index
   */
  int setCellDataArray(vtkAbstractArray *array) {
    return _geoMesh.mesh->GetCellData()->AddArray(array);
  }
  /**
   * Set field data
   * @param array array pointer
   * @return array index
   */
  int setFieldDataArray(vtkAbstractArray *array) {
    return _geoMesh.mesh->GetFieldData()->AddArray(array);
  }

  /**
   * Remove point data array by index
   * @param arrayIdx array index
   */
  void unsetPointDataArray(int arrayIdx) {
    _geoMesh.mesh->GetPointData()->RemoveArray(arrayIdx);
  }
  /**
   * Remove point data array by name
   * @param arrayName array name
   */
  void unsetPointDataArray(const std::string &arrayName) {
    _geoMesh.mesh->GetPointData()->RemoveArray(arrayName.c_str());
  }
  /**
   * Remove cell data array by index
   * @param arrayIdx array index
   */
  void unsetCellDataArray(int arrayIdx) {
    _geoMesh.mesh->GetCellData()->RemoveArray(arrayIdx);
  }
  /**
   * Remove cell data array by name
   * @param arrayName array name
   */
  void unsetCellDataArray(const std::string &arrayName) {
    _geoMesh.mesh->GetCellData()->RemoveArray(arrayName.c_str());
  }
  /**
   * Remove field data array by index
   * @param arrayIdx array index
   */
  void unsetFieldDataArray(int arrayIdx) {
    _geoMesh.mesh->GetFieldData()->RemoveArray(arrayIdx);
  }
  /**
   * Remove field data array by name
   * @param arrayName array name
   */
  void unsetFieldDataArray(const std::string &arrayName) {
    _geoMesh.mesh->GetFieldData()->RemoveArray(arrayName.c_str());
  }

 protected:
  /**
   * Get the underlying geometry object.
   * @note This is not a copy; Gmsh does not support cloning a model.
   */
  const GeoMesh &getGeoMesh() const { return _geoMesh; }

 private:
  void reconstructGeo();

 private:
  GeoMesh _geoMesh;

  /// Dihedral angle threshold (in radians) to classify surfaces (Default: 30
  /// degrees)
  double _angleThreshold;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_GEOMESHBASE_H_
