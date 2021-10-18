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
#include <vtkDataObject.h>
#include <vtkFieldData.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

namespace NEM {
namespace MSH {

using nemId_t = std::int64_t;

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
class NEMOSYS_EXPORT geoMeshBase : public vtkDataObject {
  vtkAbstractTypeMacro(geoMeshBase, vtkDataObject)

 public:
  /**
   * Query for conversion from vtk to gmsh
   * @param vtkType
   * @return type given by @c GMSH_MSH_TYPES enum
   */
  static int getGmshTypeFromVTKType(int vtkType);

 public:
  geoMeshBase();
  ~geoMeshBase() override;
  // Copy and move constructors/assignment operators deleted in vtkDataObject

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
  /**
   * Take the GeoMesh of another @c geoMeshBase. Note that @p otherGeoMesh will
   * be left with an empty mesh.
   * @param otherGeoMesh other geoMeshBase object; mesh will be left empty
   */
  virtual void takeGeoMesh(geoMeshBase *otherGeoMesh);
  /**
   * Construct the geometry from the mesh alone. Sets geometry, link, and side
   * set. Does not alter the mesh.
   */
  virtual void reconstructGeo();

  /**
   * Merge mesh data from a new geoMesh to an existing geoMesh. Deletes the data
   * from the new geoMesh pointer
   * @param otherGeoMesh A new geoMesh data to merge
   */
  virtual void mergeGeoMesh(geoMeshBase *otherGeoMesh);

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
    getPoint(id, &x);
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
   * Get cell.
   * @param cellId cell id
   * @return cell
   */
  vtkCell *getCell(nemId_t cellId) const {
    return _geoMesh.mesh->GetCell(cellId);
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
    getCellBounds(cellId, &bounds);
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
   * Get list of cells sharing points @p ptIds excluding @p cellId.
   * @param cellId id of cell to exclude
   * @param ptIds point id
   * @param cellIds list of cell ids
   */
  void getCellNeighbors(nemId_t cellId, vtkIdList *ptIds,
                        vtkIdList *cellIds) const {
    _geoMesh.mesh->GetCellNeighbors(cellId, ptIds, cellIds);
  }
  /**
   * Get number of arrays in the point data.
   * @return number of arrays in the point data
   */
  int getNumberOfPointDataArrays() const {
    return _geoMesh.mesh->GetPointData()->GetNumberOfArrays();
  }
  /**
   * Get copy of point data array by name. Optionally, get array index.
   * @param arrayName array name
   * @param array point data array
   * @param arrayIdx array index
   */
  void getPointDataArrayCopy(const std::string &arrayName,
                             vtkAbstractArray *array,
                             int *arrayIdx = nullptr) const;
  /**
   * Create and return a copy of point data array by name. Optionally, get array
   * index.
   * @param arrayName array name
   * @param arrayIdx array index
   * @return pointer to copy of array
   */
  vtkSmartPointer<vtkAbstractArray> getPointDataArrayCopy(
      const std::string &arrayName, int *arrayIdx = nullptr) const;
  /**
   * Get copy of point data array by index.
   * @param arrayIdx array index
   * @param array point data array
   */
  void getPointDataArrayCopy(int arrayIdx, vtkAbstractArray *array) const;
  /**
   * Create and return a copy of point data array by index.
   * @param arrayIdx array index
   * @return pointer to copy of array
   */
  vtkSmartPointer<vtkAbstractArray> getPointDataArrayCopy(int arrayIdx) const;
  /**
   * Get number of arrays in the cell data.
   * @return number of arrays in the cell data
   */
  int getNumberOfCellDataArrays() const {
    return _geoMesh.mesh->GetCellData()->GetNumberOfArrays();
  }
  /**
   * Get copy of cell data array by name. Optionally, get array index.
   * @param arrayName array name
   * @param array cell data array
   * @param arrayIdx array index
   */
  void getCellDataArrayCopy(const std::string &arrayName,
                            vtkAbstractArray *array,
                            int *arrayIdx = nullptr) const;
  /**
   * Create and return a copy of cell data array by name. Optionally, get array
   * index.
   * @param arrayName array name
   * @param arrayIdx array index
   * @return pointer to copy of array
   */
  vtkSmartPointer<vtkAbstractArray> getCellDataArrayCopy(
      const std::string &arrayName, int *arrayIdx = nullptr) const;
  /**
   * Get copy of cell data array by index.
   * @param arrayIdx array index
   * @param array cell data array
   */
  void getCellDataArrayCopy(int arrayIdx, vtkAbstractArray *array) const;
  /**
   * Create and return a copy of cell data array by index.
   * @param arrayIdx array index
   * @return pointer to copy of array
   */
  vtkSmartPointer<vtkAbstractArray> getCellDataArrayCopy(int arrayIdx) const;
  /**
   * Get number of arrays in the field data.
   * @return number of arrays in the field data
   */
  int getNumberOfFieldDataArrays() const {
    return _geoMesh.mesh->GetFieldData()->GetNumberOfArrays();
  }
  /**
   * Get copy of field data array by name. Optionally, get array index.
   * @param arrayName array name
   * @param array field data array
   * @param arrayIdx array index
   */
  void getFieldDataArrayCopy(const std::string &arrayName,
                             vtkAbstractArray *array,
                             int *arrayIdx = nullptr) const;
  /**
   * Create and return a copy of field data array by name. Optionally, get array
   * index.
   * @param arrayName array name
   * @param arrayIdx array index
   * @return pointer to copy of array
   */
  vtkSmartPointer<vtkAbstractArray> getFieldDataArrayCopy(
      const std::string &arrayName, int *arrayIdx = nullptr) const;
  /**
   * Get copy of field data array by index.
   * @param arrayIdx array index
   * @param array field data array
   */
  void getFieldDataArrayCopy(int arrayIdx, vtkAbstractArray *array) const;
  /**
   * Create and return a copy of field data array by index.
   * @param arrayIdx array index
   * @return pointer to copy of array
   */
  vtkSmartPointer<vtkAbstractArray> getFieldDataArrayCopy(int arrayIdx) const;

 protected:
  struct NEMOSYS_EXPORT SideSet {
    SideSet(vtkPolyData *sideSet, vtkIntArray *geoEnt,
            vtkIdTypeArray *origCell = nullptr, vtkIntArray *cellFace = nullptr,
            vtkStringArray *setNames = nullptr);
    SideSet() = default;
    explicit SideSet(vtkPolyData *sides);
    /**
     * @brief Cells represent edges/faces of some GeoMesh.
     */
    vtkSmartPointer<vtkPolyData> sides{nullptr};

    vtkSmartPointer<vtkIntArray> getGeoEntArr() const;
    void setGeoEntArr(vtkIntArray *arr);
    vtkSmartPointer<vtkIdTypeArray> getOrigCellArr() const;
    void setOrigCellArr(vtkIdTypeArray *arr);
    vtkSmartPointer<vtkIntArray> getCellFaceArr() const;
    void setCellFaceArr(vtkIntArray *arr);
    vtkSmartPointer<vtkStringArray> getSideSetNames() const;
    void setSideSetNamesArr(vtkStringArray *arr);

   private:
    static constexpr auto GEO_ENT_NAME = "GeoEnt";
    // Set to match dataSetRegionSurfaceFilter
    static constexpr auto ORIG_CELL_NAME = "OrigCellIds";
    static constexpr auto CELL_FACE_NAME = "CellFaceIds";
    static constexpr auto NAME_ARR_NAME = "Side Set Names";
  };

  struct NEMOSYS_EXPORT GeoMesh {
    vtkSmartPointer<vtkUnstructuredGrid> mesh;
    // name of gmsh model with geo
    std::string geo;
    // name of cell data array in mesh (and sideSet if sideSet not null) storing
    // physical group (matching geo) for each cell
    std::string link;
    // boundaries of regions, represented by lower dimensional cells
    SideSet sideSet;

    /**
     * @brief For each edge/polygon in sideSet, find the corresponding mesh
     * element in mesh.
     * @details Does nothing if "OrigCellIds"/"CellFaceIds" arrays already
     * present. Sets the sideSet's "OrigCellIds"/"CellFaceIds" arrays.
     * Orientation of side only checked if multiple candidate original cells.
     * Both "OrigCellIds" and "CellFaceIds" are set to -1 for any cells in the
     * sideSet that can't be matched.
     */
    void findSide2OrigCell();
    /**
     * @brief Get dimension of mesh.
     * @details All cells in mesh assumed to have same dimension.
     * @return dimension of mesh
     */
    int getDimension() const;
  };

  static constexpr auto GEO_ENT_DEFAULT_NAME = "GeoEnt";

 protected:
  /**
   * Construct geoMeshBase from a GeoMesh struct
   * @param inGeoMesh a mesh, geo, link, and sideSet
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
  /**
   * Set the underlying geometry object.
   * @param geoMesh a mesh and geometry
   */
  void setGeoMesh(const geoMeshBase::GeoMesh &geoMesh) {
    assert(geoMesh.mesh->GetNumberOfCells() == 0 ||
           geoMesh.mesh->GetNumberOfPoints() == 0 ||
           geoMesh.sideSet.sides == nullptr ||
           geoMesh.sideSet.sides->GetPoints() == geoMesh.mesh->GetPoints());
    _geoMesh = geoMesh;
  }

 private:
  virtual void resetNative() = 0;

 private:
  GeoMesh _geoMesh;

  /// Dihedral angle threshold (in radians) to classify surfaces (Default: 30
  /// degrees)
  double _angleThreshold;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_GEOMESHBASE_H_
