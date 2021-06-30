#ifndef meshBase_H
#define meshBase_H

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <vtkCellLocator.h>
#include <vtkDataSet.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkModelMetadata.h>
#include <vtkSmartPointer.h>
#include <vtkStaticCellLocator.h>
#include <vtkStaticPointLocator.h>
#include <vtkUnstructuredGrid.h>

#include "nemosys_export.h"

typedef std::size_t nemId_t;

/** @brief A brief description of meshBase.
    @note Virtual methods are usually implemented in vtkMesh.C. We use that
        class for more VTK-specific functions that we want wrapped by meshBase
    @note Try using smart pointer instances of meshBase in drivers and pass the
        underlying raw pointer with the get member function of the smart
        pointer when required (e.g., source->transfer(target.get(), ...))
        for smart pointer instance of target. Ideally, driver classes should
        have empty destructors, which can be realized by using smart pointers.
        In general, smart pointer should be used and scoped when they can.
**/

class NEMOSYS_DEPRECATED_EXPORT meshBase {
  // --- constructors and destructors
 public:
  meshBase()
      : numPoints(0),
        numCells(0),
        dataSet(nullptr),
        checkQuality(false),
        continuous(false),
        metadata(nullptr) {
    std::cout << "meshBase constructed" << std::endl;
  }

  //    meshBase(const meshBase &mb) = delete;
  //    meshBase &operator=(const meshBase &mb) = delete;
  //    meshBase(meshBase &&mb) = delete;
  //    meshBase &operator=(meshBase &&mb) = delete;

  // Dropping support for copy constructor
  // meshBase(const meshBase& mb)
  //{
  //    numPoints = mb.numPoints;
  //    numCells = mb.numCells;
  //    hasSizeField = mb.hasSizeField;
  //    dataSet = mb.dataSet;
  //    filename = mb.filename;
  //    checkQuality = mb.checkQuality;
  //    order = mb.order;
  //    newArrayNames = mb.newArrayNames;
  //}

  virtual ~meshBase() { std::cout << "meshBase destroyed" << std::endl; }

  // --- meshBase factory methods
 public:
  /** @brief Construct vtkMesh from filename.
      @param fname name of mesh file
      @return <>
  **/
  static meshBase *Create(const std::string &fname);

  /** @brief Construct from existing vtkDataSet and assign newname as
          filename.
      @param other The vtkDataSet used to construct the mesh
      @param fname name of mesh file
      @return <>
  **/
  static meshBase *Create(vtkSmartPointer<vtkDataSet> other,
                          const std::string &newname);

  /** @brief create from coordinates and connectivities.
      @param xCrds <>
      @param yCrds <>
      @param zCrds <>
      @param elmConn <>
      @param cellType one of the vtkCellType enums. Currently, only VTK_TETRA
          and VTK_TRIANGLE are supported.
      @param newname name of mesh file
      @return <>
  **/
  static meshBase *Create(const std::vector<double> &xCrds,
                          const std::vector<double> &yCrds,
                          const std::vector<double> &zCrds,
                          const std::vector<nemId_t> &elmConn,
                          const int cellType, const std::string &newname);

  /** @brief Create shared ptr from fname.
      @param fname name of mesh file
      @return <>
  **/
  static std::shared_ptr<meshBase> CreateShared(const std::string &fname);

  /** @brief Create shared ptr from existing meshbase.
      @param mesh the existing meshbase
      @return <>
  **/
  static std::shared_ptr<meshBase> CreateShared(meshBase *mesh);

  /** @brief Create shared ptr from existing vtkDataset and assign newname as
          filename.
      @param other The vtkDataSet used to construct the mesh
      @param fname name of mesh file
      @return <>
  **/
  static std::shared_ptr<meshBase> CreateShared(
      vtkSmartPointer<vtkDataSet> other, const std::string &newname);

  /** @brief Version of raw data mesh creation for memory managed shared_ptr
          instance.
      @param xCrds <>
      @param yCrds <>
      @param zCrds <>
      @param elmConn <>
      @param cellType one of the vtkCellType enums. Currently, only
          VTK_TETRA and VTK_TRIANGLE are supported.
      @param fname name of mesh file
      @return <>
  **/
  static std::shared_ptr<meshBase> CreateShared(
      const std::vector<double> &xCrds, const std::vector<double> &yCrds,
      const std::vector<double> &zCrds, const std::vector<nemId_t> &elmConn,
      int cellType, const std::string &newname);

  /** @brief create unique ptr from fname
      @param fname name of mesh file
      @return <>
  **/
  static std::unique_ptr<meshBase> CreateUnique(const std::string &fname);

  /** @brief version of raw data mesh creation for memory managed unique ptr
          instance
      @param xCrds <>
      @param yCrds <>
      @param zCrds <>
      @param elmConn <>
      @param cellType one of the vtkCellType enums. Currently, only VTK_TETRA
          and VTK_TRIANGLE are supported.
      @param fname name of mesh file
      @return <>
  **/
  static std::unique_ptr<meshBase> CreateUnique(
      const std::vector<double> &xCrds, const std::vector<double> &yCrds,
      const std::vector<double> &zCrds, const std::vector<nemId_t> &elmConn,
      int cellType, const std::string &newname);
  /** @brief construct from existing vtkDataSet and assign newname as filename
      @param other The vtkDataSet used to construct the mesh
      @param fname name of mesh file
      @return <>
  **/
  static std::unique_ptr<meshBase> CreateUnique(
      vtkSmartPointer<vtkDataSet> other, const std::string &newname);
  /** @brief construct from existing meshbase object
      @param mesh the existing meshbase
      @return <>
  **/
  static std::unique_ptr<meshBase> CreateUnique(meshBase *mesh);

  /** @brief construct vtkMesh from gmsh msh file (called in Create methods)
      @param fname name of mesh file
      @return <>
  **/
  static meshBase *exportGmshToVtk(const std::string &fname);

  /** @brief construct vtkMesh from netgen vol file (called in Create methods)
      @param fname name of mesh file
      @return <>
  **/
  static meshBase *exportVolToVtk(const std::string &fname);

  /** @brief construct vtkMesh from netgen vol file (called in Create methods)
      @param fname name of mesh file
      @return <>
  **/
  static meshBase *exportPntToVtk(const std::string &fname);

  /** @brief construct vtkMesh from exodusII files
      @param fname name of mesh file
      @return <>
  **/
  static meshBase *exportExoToVtk(const std::string &fname);

  /*
  static meshBase *generateMesh(const std::string &fname,
                                const std::string &meshEngine,
                                meshingParams *params);
  */

  /** @brief stitch together several meshBases
      @param mbObjs a vector of meshBase objects to stich together
      @return <>
  **/
  static meshBase *stitchMB(const std::vector<meshBase *> &mbObjs);

  /** @brief stitch together several meshBase
      @param mbObjs a vector of meshBase objects to stich together
      @return <>
  **/
  static std::shared_ptr<meshBase> stitchMB(
      const std::vector<std::shared_ptr<meshBase>> &_mbObjs);

  /** @brief mesh partitioning (with METIS)
      @param mbObj The meshBase object to partition.
      @param numPartitions The number of partitions to partition the mesh into
      @return <>
  **/
  static std::vector<std::shared_ptr<meshBase>> partition(const meshBase *mbObj,
                                                          int numPartitions);

  /** @brief extract subset of mesh given list of cell ids and return meshBase
          obj
      @param mesh The meshBase object to extract the subset from.
      @param cellIds <>
      @return meshBase object representing the subset
  **/
  static meshBase *extractSelectedCells(meshBase *mesh,
                                        const std::vector<nemId_t> &cellIds);

  /** @brief helper wrapped by function above
      @param mesh The meshBase object to extract the subset from.
      @param cellIds <>
      @return meshBase object representing the subset
  **/
  static meshBase *extractSelectedCells(
      vtkSmartPointer<vtkDataSet> mesh,
      vtkSmartPointer<vtkIdTypeArray> cellIds);

  // --- access
 public:
  /** @brief abstract read method reserved for derived classes
      @param fname name of mesh file
  **/
  virtual void read(const std::string &fname) = 0;

  /** @brief get point with id
      @param id The id of the point.
      @return <>
  **/
  virtual std::vector<double> getPoint(nemId_t id) const = 0;

  /** @brief get 3 vecs with x,y and z coords
      @return <>
  **/
  virtual std::vector<std::vector<double>> getVertCrds() const = 0;

  /** @brief get cell with id
      @param id The id of the cell.
      @return point indices and respective coordinates
  **/
  virtual std::map<nemId_t, std::vector<double>> getCell(nemId_t id) const = 0;

  /** @brief get vector of coords of cell with id
      @param id The id of the cell.
      @return vector of coords of cell
  **/
  virtual std::vector<std::vector<double>> getCellVec(nemId_t id) const = 0;

  /** @brief get edge lengths of dataSet
      @param ofname <>
  **/
  virtual void inspectEdges(const std::string &ofname) const = 0;

  /** @brief get this meshes' dataSet
      @return this mesh's dataSet
  **/
  vtkSmartPointer<vtkDataSet> getDataSet() const { return dataSet; }

  /** @brief extract the surface mesh
      @return the surface mesh for this mesh.
  **/
  virtual vtkSmartPointer<vtkDataSet> extractSurface() = 0;

  /** @brief register data to dataSet's point data
      @param name <>
      @param data <>
  **/
  virtual void setPointDataArray(const std::string &name,
                                 const std::vector<std::vector<double>> &data) {
  }

  /** @brief register data to dataSet's point data
      @param name <>
      @param data <>
 **/
  virtual void setPointDataArray(const std::string &name,
                                 const std::vector<double> &data) {}

  /** @brief register data to dataSet's cell data
      @param name <>
      @param data <>
  **/
  virtual void setCellDataArray(const std::string &name,
                                const std::vector<std::vector<double>> &data) {}

  /** @brief register data to dataSet's cell data
      @param name <>
      @param data <>
  **/
  virtual void setCellDataArray(const std::string &name,
                                const std::vector<double> &data) {}

  /** @brief get *scalar* point or cell data array.
          assumes data is not allocated prior to calling
      @param name <>
      @param data <>
  **/
  virtual void getPointDataArray(const std::string &name,
                                 std::vector<double> &data) {}

  /** @brief get *scalar* point or cell data array.
          assumes data is not allocated prior to calling
      @param arrayId <>
      @param data <>
  **/
  virtual void getPointDataArray(int arrayId, std::vector<double> &data) {}

  /** @brief <>
      @param name <>
      @return <>
  **/
  virtual int getCellDataIdx(const std::string &name) { return 0; }

  /** @brief <>
      @param name <>
      @param data <>
  **/
  virtual void getCellDataArray(const std::string &name,
                                std::vector<double> &data) {}

  /** @brief <>
      @param arrayId <>
      @param data <>
  **/
  virtual void getCellDataArray(int arrayId, std::vector<double> &data) {}

  /** @brief delete array with id from dataSet's point data
      @param arrayID <>
  **/
  virtual void unsetPointDataArray(int arrayID) {}

  /** @brief <>
      @param name <>
  **/
  virtual void unsetPointDataArray(const std::string &name) {}

  /** @brief delete array with id from dataSet's cell data
      @param arrayID <>
  **/
  virtual void unsetCellDataArray(int arrayID) {}

  /** @brief <>
      @param name <>
  **/
  virtual void unsetCellDataArray(const std::string &name) {}

  /** @brief delete array with id from dataSet's field data
      @param name <>
  **/
  virtual void unsetFieldDataArray(const std::string &name) {}

  /** @brief get diameter of circumsphere of each cell
      @return <>
  **/
  virtual std::vector<double> getCellLengths() const = 0;

  /** @brief get center of a cell
      @param cellID <>
      @return <>
  **/
  virtual std::vector<double> getCellCenter(nemId_t cellID) const = 0;

  /** @brief build locators for efficient search operations
      @return <>
  **/
  vtkSmartPointer<vtkStaticCellLocator> buildStaticCellLocator();

  /** @brief build thread-safe point locator for efficient search operations
   *  @return <>
   */
  vtkSmartPointer<vtkStaticPointLocator> buildStaticPointLocator();

  /** @brief get cell type as an integer
      assumes all elements are the same type
      @return <>
  **/
  virtual int getCellType() const = 0;

  /** @brief get connectivities.

      This is only safe to use if mesh has cells of the same type or you have
      information on the number of cells of each type and
      the order in which they appear (for look up in resulting vector)

      @return <>
  **/
  virtual std::vector<nemId_t> getConnectivities() const = 0;
  // set metadata, including sidesets
  void setMetadata(vtkSmartPointer<vtkModelMetadata> _metadata) {
    metadata = _metadata;
  }
  vtkSmartPointer<vtkModelMetadata> getMetadata() { return metadata; }

  // --- integration
 public:
  /** @brief integrate arrays in arrayIDs over the mesh.
      @param arrayIDs <>
      @return total integral for each datum is returned
  **/
  std::vector<std::vector<double>> integrateOverMesh(
      const std::vector<int> &arrayIDs);

  // --- size field generation
 public:
  /** @brief generate size field based on method and given a point data array.
      @param method (e.g., "gradient", "value", "error estimator")
      @param arrayID <>
      @param dev_mlt used to determine which cells to consider for refinement
      @param maxIsmin used to determine which cells to consider for refinement
      @param sizeFactor <>
      @param order <>
  **/
  void generateSizeField(const std::string &method, int arrayID, double dev_mlt,
                         bool maxIsmin, double sizeFactor = 1.0, int order = 1);

  /** @brief check for named array in vtk and return its integer id.
      @param pointOrCell boolean that tells the method whether to transfer
          point (False) or cell (True) data.
      @return <>
  **/
  int IsArrayName(const std::string &name, bool pointOrCell = false) const;

  // --- adaptive mesh refinement
 public:
  /** @brief perform sizefield-based h-refinement.
      @param method <>
      @param arrayID <>
      @param dev_mult <>
      @param maxIsmin <>
      @param edge_scale <>
      @param ofname <>
      @param transferData <>
      @param sizeFactor <>
  **/
  void refineMesh(const std::string &method, int arrayID, double dev_mult,
                  bool maxIsmin, double edge_scale, const std::string &ofname,
                  bool transferData, double sizeFactor = 1.,
                  bool constrainBoundary = false);

  /** @brief perform sizefield-based h-refinement.
      @param method <>
      @param arrayName <>
      @param dev_mult <>
      @param maxIsmin <>
      @param edge_scale <>
      @param ofname <>
      @param transferData <>
      @param sizeFactor <>
  **/
  void refineMesh(const std::string &method, const std::string &arrayName,
                  double dev_mult, bool maxIsmin, double edge_scale,
                  const std::string &ofname, bool transferData,
                  double sizeFactor = 1.);

  /** @brief added for uniform refinement by driver
      @param method <>
      @param edge_scale <>
      @param ofname <>
      @param transferData <>
  **/
  void refineMesh(const std::string &method, double edge_scale,
                  const std::string &ofname, bool transferData,
                  bool constrainBoundary = false);

  /** @brief <>
      @param method <>
      @param arrayID <>
      @param order <>
      @param ofname <>
      @param transferData <>
  **/
  void refineMesh(const std::string &method, int arrayID, int order,
                  const std::string &ofname, bool transferData);

  /** @brief <>
      @param method <>
      @param arrayName <>
      @param order <>
      @param ofname <>
      @param transferData <>
  **/
  void refineMesh(const std::string &method, const std::string &arrayName,
                  int order, const std::string &ofname, bool transferData);

  // --- diagnostics
 public:
  /** @brief generate a report of the mesh
   **/
  virtual void report() const {};

  /** @brief return the number of points
      @return The number of points.
  **/
  nemId_t getNumberOfPoints() const { return numPoints; }

  /** @brief return the number of cells
      @return The number of cells.
  **/
  nemId_t getNumberOfCells() const { return numCells; }

  /** @brief <>
      @param ofname <>
  **/
  void checkMesh(const std::string &ofname) const;

  // --- for distributed data sets.
 public:
  /** @brief global to local mapping of nodes
      @note These are only generated if mesh is one of the partitions returned
           from a call to meshBase::partition
  **/
  std::map<nemId_t, nemId_t> getGlobToPartNodeMap() {
    return globToPartNodeMap;
  }

  /** @brief global to local mapping of cells
      @note These are only generated if mesh is one of the partitions returned
           from a call to meshBase::partition
  **/
  std::map<nemId_t, nemId_t> getGlobToPartCellMap() {
    return globToPartCellMap;
  }

  /** @brief local to global mapping of nodes
      @note These are only generated if mesh is one of the partitions returned
           from a call to meshBase::partition
  **/
  std::map<nemId_t, nemId_t> getPartToGlobNodeMap() {
    return partToGlobNodeMap;
  }

  /** @brief local to global mapping of cells
      @note These are only generated if mesh is one of the partitions returned
           from a call to meshBase::partition
  **/
  std::map<nemId_t, nemId_t> getPartToGlobCellMap() {
    return partToGlobCellMap;
  }

  // --- write and conversion
 public:
  /** @brief write the mesh to file named after the private var 'filename'.

      The file extension of the private var "filename" determines the format of
  the output file
  **/
  virtual void write() const { write(filename); }

  /** @brief write the mesh to file named fname
      @param fname The name of the file to write to
  **/
  virtual void write(const std::string &fname) const = 0;

  /** @brief convert to gmsh format without data
      @param outputStream <>
  **/
  void writeMSH(std::ofstream &outputStream);

  /** @brief convert to gmsh format without data
      @param fname The name of the file to write to
  **/
  void writeMSH(const std::string &fname);

  /** @brief convert to gmsh format with specified point or cell data
      @param outputStream <>
      @param pointOrCell <>
      @param arrayID <>
  **/
  void writeMSH(std::ofstream &outputStream, const std::string &pointOrCell,
                int arrayID);

  /** @brief convert to gmsh format without data
      @param fname The name of the file to write to
      @param pointOrCell <>
      @param arrayID <>
  **/
  void writeMSH(const std::string &fname, const std::string &pointOrCell,
                int arrayID);

  /** @brief convert to gmsh format with specified point or cell data for
          only volume elements (USE ONLY FOR MADLIB STUFF)
      @param outputStream <>
      @param pointOrCell <>
      @param arrayID <>
      @param onlyVol <>
  **/
  void writeMSH(std::ofstream &outputStream, const std::string &pointOrCell,
                int arrayID,
                bool onlyVol);  // added for overloading, doesn't do anything

  /** @brief convert to gmsh format with specified point or cell data for
          only volume elements (USE ONLY FOR MADLIB STUFF)
      @param fname The name of the file to write to
      @param pointOrCell <>
      @param arrayID <>
      @param onlyVol <>
  **/
  void writeMSH(const std::string &fname, const std::string &pointOrCell,
                int arrayID, bool onlyVol);

  /** @brief surfWithPatch must have patchNo array
      @param surfWithPatch <>
      @param mapFile <>
      @param outputStream <>
  **/
  void writeCobalt(meshBase *surfWithPatch, const std::string &mapFile,
                   std::ofstream &outputStream);

  /** @brief surfWithPatch must have patchNo array
      @param surfWithPatch <>
      @param mapFile <>
      @param ofname <>
  **/
  void writeCobalt(meshBase *surfWithPatch, const std::string &mapFile,
                   const std::string &ofname);

  /** @brief set the file name.

      This will allow vtk to dispatch appropriate writers
      based on the extension and whether it is supported by vtk.

      @param fname The name to set the private variable "filename" to
  **/
  void setFileName(const std::string &fname) { filename = fname; }

  /** @brief get the current file name
      @return The current value of the private variable "filename"
  **/
  const std::string &getFileName() const { return filename; }

  /** @brief set whether to check quality of transfer by back-transfer and rmse
      @param x <>
  **/
  void setCheckQuality(bool x) { checkQuality = x; }

  /** @brief set weighted averaging/smoothing for cell data transfer
          (default is off)
      @param x <>
  **/
  void setContBool(bool x) { continuous = x; }

  // remove all quad elements by naive conversion to tris
  meshBase *convertQuads();

  /** @brief get new array names for use in transfer
   **/
  std::vector<std::string> getNewArrayNames() { return newArrayNames; }

  /** @brief given array names, return corresponding ids
   **/
  std::vector<int> getArrayIDs(std::vector<std::string> arrayNames,
                               bool fromPointArrays = false);

  /** @brief Converts given hexahedral VTK dataset into tetrahedral mesh
              and stores it into dataSet variable.
      @param meshdataSet Input hexahedral mesh dataset
  **/
  void convertHexToTetVTK(vtkSmartPointer<vtkDataSet> meshdataSet);

 protected:
  /** @brief number of points in mesh
   **/
  nemId_t numPoints;

  /** @brief number of cells in mesh
   **/
  nemId_t numCells;

  /** @brief mesh points, topology and data
   **/
  vtkSmartPointer<vtkDataSet> dataSet;

  /** @brief name of mesh file
   **/
  std::string filename;

  /** @brief check transfer quality when on
   **/
  bool checkQuality;

  /** @brief switch on / off weighted averaging for cell data transfer
          (default is off)
  **/
  bool continuous;

  /** @brief new names to set for transferred data
   **/
  std::vector<std::string> newArrayNames;

  /** @brief map between global and local node idx in partition
          for distributed data sets

      Only populated for mesh resulting from call to meshBase::partition
  **/
  std::map<nemId_t, nemId_t> globToPartNodeMap;

  /** @brief map between global and local cell idx in partition
   **/
  std::map<nemId_t, nemId_t> globToPartCellMap;

  /** @brief map between local and global node idx in partition
   **/
  std::map<nemId_t, nemId_t> partToGlobNodeMap;

  /** @brief map between local and global cell idx in partition
   **/
  std::map<nemId_t, nemId_t> partToGlobCellMap;
  // metadata
  vtkSmartPointer<vtkModelMetadata> metadata;
};

/** @brief sum comparison for vectors representing faces inserted into map
 **/
struct NEMOSYS_EXPORT sortNemId_tVec_compare {
  // TODO: This is a comparison operation making copies since passed by value!
  //       Should be pass by const-ref. If sorting is necessary, do it
  //       beforehand.
  bool operator()(std::vector<nemId_t> lhs, std::vector<nemId_t> rhs) const;
};

// --- auxiliary helpers
/** @brief compares two meshBase classes. used in testing
    @param mesh1 <>
    @param mesh2 <>
**/
NEMOSYS_EXPORT int diffMesh(meshBase *mesh1, meshBase *mesh2);

#ifndef SWIG
/** @brief write patch map file for roc prep (trivial identity mapping)
    @param mapFile <>
    @param patchMap <>
**/
NEMOSYS_EXPORT void writePatchMap(const std::string &mapFile,
                                  const std::map<int, int> &patchMap);

/** @brief write patch map file for roc prep (trivial identity mapping)
    @param outputStream <>
    @param patchMap <>
**/
NEMOSYS_EXPORT void writePatchMap(std::ofstream &outputStream,
                                  const std::map<int, int> &patchMap);
#endif  // SWIG

#endif  // meshBase_H
