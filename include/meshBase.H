#ifndef meshBase_H
#define meshBase_H

// VTK
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>

// Nemosys
#include <pntMesh.H>
#include <cobalt.H>
#include <patran.H>

// stl
#include <vector>
#include <iosfwd>
#include <fstream>
#include <map>
#include <string>
#include <memory>


/* NOTE: Virtual methods are usually implemented in vtkMesh.C. We use that class
         for more VTK-specific functions that we want wrapped by meshBase */
/* NOTE: Try using smart pointer instances of meshBase in drivers and pass the
         underlying raw pointer with the get member function of the smart
         pointer when required (e.g., source->transfer(target.get(), ...))
         for smart pointer instance of target. Ideally, driver classes should
         have empty destructors, which can be realized by using smart pointers.
         In general, smart pointer should be used and scoped when they can. */

class meshingParams;

class meshBase 
{
  // --- constructors and destructors
  public:
    meshBase()
        : dataSet(nullptr), numPoints(0), numCells(0),
          hasSizeField(false), checkQuality(false), continuous(false), order(1)
    {
      std::cout << "meshBase constructed" << std::endl;
    }

    // Dropping support for copy constructor
    //meshBase(const meshBase& mb)
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

    virtual ~meshBase()
    {
      std::cout << "meshBase destroyed" << std::endl;
    }

  // --- meshBase factory methods
  public:
    // construct vtkMesh from filename
    // this method calls the other factory methods based on extension
    // caller must delete object after use
    static meshBase* Create(const std::string &fname);
    // construct from existing vtkDataSet and assign newname as filename
    // caller must delete object after use
    static meshBase* Create(vtkSmartPointer<vtkDataSet> other,
                            const std::string &newname);
    /* create from coordinates and connectivities.
       NOTE: use of this is only valid when mesh has one cell type. 
             cellType one of the vtkCellType enums.
             Currently, only VTK_TETRA and VTK_TRIANGLE are supported
       Caller must delete object after use  */
    static meshBase* Create(const std::vector<double>& xCrds,
                            const std::vector<double>& yCrds,
                            const std::vector<double>& zCrds,
                            const std::vector<int>& elmConn,
                            const int cellType,
                            const std::string &newname);
    // create shared ptr from fname
    // memory is managed by shared pointer, so do not call delete after use
    static std::shared_ptr<meshBase> CreateShared(const std::string &fname);
    // create shared ptr from existing meshbase
    // memory is managed by shared pointer, so do not call delete after use
    // (be careful with this one!)
    static std::shared_ptr<meshBase> CreateShared(meshBase* mesh);
    // create shared ptr from existing vtkDataset and assign newname as filename
    // memory is managed by shared pointer, so do not call delete after use
    static std::shared_ptr<meshBase>
    CreateShared(vtkSmartPointer<vtkDataSet> other,
                 const std::string &newname);
    // version of raw data mesh creation for memory managed shared_ptr instance (do not call delete)
    static std::shared_ptr<meshBase>
    CreateShared(const std::vector<double> &xCrds,
                 const std::vector<double> &yCrds,
                 const std::vector<double> &zCrds,
                 const std::vector<int> &elmConn,
                 const int cellType,
                 const std::string &newname);
    // create unique ptr from fname
    // memory is managed by unique pointer, so do not call delete after use
    static std::unique_ptr<meshBase> CreateUnique(const std::string &fname);
    // version of raw data mesh creation for memory managed unique ptr instance (do not call delete)
    static std::unique_ptr<meshBase>
    CreateUnique(const std::vector<double> &xCrds,
                 const std::vector<double> &yCrds,
                 const std::vector<double> &zCrds,
                 const std::vector<int> &elmConn,
                 const int cellType,
                 const std::string &newname);
    // construct from existing vtkDataSet and assign newname as filename
    // memory is managed by unique ptr, so do not call delete
    static std::unique_ptr<meshBase>
    CreateUnique(vtkSmartPointer<vtkDataSet> other,
                 const std::string &newname);
    // construct from existing meshbase object
    // memory is managed by unique ptr, so do not call delete
    static std::unique_ptr<meshBase> CreateUnique(meshBase *mesh);
    // construct vtkMesh from gmsh msh file (called in Create methods)
    static meshBase *exportGmshToVtk(const std::string &fname);
    // construct vtkMesh from netgen vol file (called in Create methods)
    static meshBase *exportVolToVtk(const std::string &fname);
    // construct vtkMesh from netgen vol file (called in Create methods)
    static meshBase *exportPntToVtk(const std::string &fname);
    // construct vtkMesh from exodusII files
    static meshBase *exportExoToVtk(const std::string &fname);
    // generate volume mesh from stl
    // caller must delete object after use
    static meshBase *generateMesh(const std::string &fname,
                                  const std::string &meshEngine,
                                  meshingParams *params);
    // stitch together several meshBases
    // caller must delete object after use
    static meshBase *stitchMB(const std::vector<meshBase *> &mbObjs);
    // stitch together several meshBase
    // memory is managed by shared ptr, so do not call delete
    static std::shared_ptr<meshBase>
    stitchMB(const std::vector<std::shared_ptr<meshBase>> &_mbObjs);
    // mesh partitioning (with METIS)
    // memory is managed by shared ptr, so do not call delete
    static std::vector<std::shared_ptr<meshBase>>
    partition(const meshBase *mbObj, const int numPartitions);
    // extract subset of mesh given list of cell ids and return meshBase obj
    static meshBase *
    extractSelectedCells(meshBase *mesh,
                         const std::vector<vtkIdType> &cellIds);
    // helper wrapped by function above
    static meshBase *
    extractSelectedCells(vtkSmartPointer<vtkDataSet> mesh,
                         vtkSmartPointer<vtkIdTypeArray> cellIds);

  // --- access
  public:
    // get point with id
    virtual std::vector<double> getPoint(int id) const {}
    // get 3 vecs with x,y and z coords
    virtual std::vector<std::vector<double>> getVertCrds() const {}
    // get cell with id : returns point indices and respective coordinates
    virtual std::map<int, std::vector<double>> getCell(int id) const {}
    // get vector of coords of cell with id
    virtual std::vector<std::vector<double>> getCellVec(int id) const {}
    // get edge lengths of dataSet
    virtual void inspectEdges(const std::string &ofname) const = 0;
    // get this meshes' dataSet
    vtkSmartPointer<vtkDataSet> getDataSet() const { return dataSet; }
    // extract the surface mesh
    virtual vtkSmartPointer<vtkDataSet> extractSurface(){}
    // register data to dataSet's point data
    virtual void
    setPointDataArray(const char *name,
                      const std::vector<std::vector<double>> &data) {}
    virtual void setPointDataArray(const char *name,
                                   const std::vector<double> &data) {}
    // register data to dataSet's cell data
    virtual void
    setCellDataArray(const char *name,
                     const std::vector<std::vector<double>> &data) {}
    virtual void setCellDataArray(const char *name,
                                  const std::vector<double> &data) {}
    // get *scalar* point or cell data array.
    // assumes data is not allocated prior to calling
    virtual void getPointDataArray(const std::string &name,
                                   std::vector<double> &data) {}
    virtual void getPointDataArray(int arrayId,
                                   std::vector<double> &data) {}
    virtual int getCellDataIdx(const char *name) {}
    virtual void getCellDataArray(const std::string &name,
                                  std::vector<double> &data) {}
    virtual void getCellDataArray(int arrayId,
                                  std::vector<double> &data) {}
    // delete array with id from dataSet's point data
    virtual void unsetPointDataArray(int arrayID) {}
    virtual void unsetPointDataArray(const char *name) {}
    // delete array with id from dataSet's cell data
    virtual void unsetCellDataArray(int arrayID) {}
    virtual void unsetCellDataArray(const char *name) {}
    // delete array with id from dataSet's field data
    virtual void unsetFieldDataArray(const char *name) {}
    // get diameter of circumsphere of each cell
    virtual std::vector<double> getCellLengths() const {}
    // get center of a cell
    virtual std::vector<double> getCellCenter(int cellID) const {}
    // build locators for efficient search operations
    vtkSmartPointer<vtkCellLocator> buildLocator();
    // get cell type as an integer
    // assumes all elements are the same type
    virtual int getCellType() const = 0;
    // get connectivities. this is only safe to use if mesh has cells of the same type
    // or you have information on the number of cells of each type and the order in 
    // which they appear (for look up in resulting vector)
    virtual std::vector<int> getConnectivities() const {}

  // --- data transfer
  public:
    /* transfer point data (pointOrCell=0) or cell data (pointOrCell=1)
       with given id from this mesh to target
       method can be "Consistent Interpolation", "Mortar Element", "RBF", etc.
       only "Consistent Interpolation" has been implemented */
    int transfer(meshBase *target, const std::string &method,
                 const std::vector<int> &arrayIDs, bool pointOrCell = false);
    int transfer(meshBase *target,
                 const std::string &method,
                 const std::vector<std::string> &arrayNames,
                 bool pointOrCell = false);
    // transfer all point and cell data from this mesh to target
    int transfer(meshBase *target, const std::string &method);

  // --- integration
  public:
    // integrate arrays in arrayIDs over the mesh. 
    // integrated data is available per cell as vtk cell data after operation
    // total integral for each datum is returned
    std::vector<std::vector<double>>
    integrateOverMesh(const std::vector<int> &arrayIDs);
 
  // --- size field generation
  public:
    // generate size field based on method (e.g., "gradient", "value",
    // "error estimator") and given a point data array.
    // dev_mlt and maxIsmin are used to determine which cells
    // to consider for refinement
    void generateSizeField(const std::string &method, int arrayID,
                           double dev_mlt, bool maxIsmin,
                           double sizeFactor = 1.);
    // set whether this mesh has a sizefield
    void setSFBool(bool q) { hasSizeField = q; }
    // get this meshes' sizefield
    bool getSFBool() const { return hasSizeField; }
    // check for named array in vtk and return its integer id (0-point data, 1-cell data)
    int IsArrayName(const std::string &name, bool pointOrCell = false) const;
    // set element shape function order (default is 1 when meshBase is constructed)
    void setOrder(int _order) { order = _order; }
    // get the shape function order
    int getOrder() const { return order; }
 
  // --- adaptive mesh refinement
  public:
    // perform sizefield-based h-refinement. edge_scale is for uniform
    // refinement and is ignored in calls where method is "gradient", "value",
    // etc. instead of "uniform"
    void refineMesh(const std::string &method, int arrayID,
                    double dev_mult, bool maxIsmin, double edge_scale,
                    const std::string &ofname, bool transferData,
                    double sizeFactor = 1.);
    void refineMesh(const std::string &method, const std::string &arrayName,
                    double dev_mult, bool maxIsmin, double edge_scale,
                    const std::string &ofname, bool transferData,
                    double sizeFactor = 1.);
    // added for uniform refinement by driver
    void refineMesh(const std::string &method, double edge_scale,
                    const std::string &ofname, bool transferData);
    void refineMesh(const std::string &method,
                    int arrayID, int order,
                    const std::string &ofname, bool transferData);
    void refineMesh(const std::string &method,
                    const std::string &arrayName, int order,
                    const std::string &ofname, bool transferData);

  // --- diagnostics
  public:
    // generate a report of the mesh
    virtual void report() const {};
    // return the number of points
    vtkIdType getNumberOfPoints() const { return numPoints; }
    // return the number of cells
    vtkIdType getNumberOfCells() const { return numCells; }
    void checkMesh(const std::string &ofname) const;

  // --- for distributed data sets. 
  public:
    /* NOTE: These are only generated if mesh is one of the partitions returned
             from a call to meshBase::partition */
    // global to local mapping of nodes
    std::map<vtkIdType, vtkIdType> getGlobToPartNodeMap() {
      return globToPartNodeMap;
    }
    // global to local mapping of cells
    std::map<vtkIdType, vtkIdType> getGlobToPartCellMap() {
      return globToPartCellMap;
    }
    // local to global mapping of nodes
    std::map<vtkIdType, vtkIdType> getPartToGlobNodeMap() {
      return partToGlobNodeMap;
    }
    // local to global mapping of cells
    std::map<vtkIdType, vtkIdType> getPartToGlobCellMap() {
      return partToGlobCellMap;
    }

  // --- write and conversion
  public:
    // write the mesh to file named after the private var 'filename'. the file
    // extension of the private var "filename" determines the format of the
    // output file
    virtual void write() const { write(filename); }
    // write the mesh to file named fname
    virtual void write(const std::string &fname) const = 0;
    // convert to gmsh format without data 
    void writeMSH(std::ofstream &outputStream);
    void writeMSH(const std::string &fname);
    // convert to gmsh format with specified point or cell data
    void writeMSH(std::ofstream &outputStream,
                  const std::string &pointOrCell, int arrayID);
    void writeMSH(const std::string &fname,
                  const std::string &pointOrCell, int arrayID);
    // convert to gmsh format with specified point or cell data for
    // only volume elements (USE ONLY FOR MADLIB STUFF)
    void writeMSH(std::ofstream &outputStream,
                  const std::string &pointOrCell, int arrayID,
                  bool onlyVol); // added for overloading, doesn't do anything 
    void writeMSH(const std::string &fname,
                  const std::string &pointOrCell, int arrayID,
                  bool onlyVol);
    // for rocstar restart hack through rflupart/prep
    // surfWithPatch must have patchNo array
    void writeCobalt(meshBase *surfWithPatch, const std::string &mapFile,
                     std::ofstream &outputStream);
    void writeCobalt(meshBase *surfWithPatch,
                     const std::string &mapFile, const std::string &ofname);
    // set the file name - this will allow vtk to dispatch appropriate writers
    // based on the extension and whether it is supported by vtk.
    void setFileName(const std::string &fname) { filename = fname; }
    // get the current file name
    const std::string &getFileName() const { return filename; }
    // set whether to check quality of transfer by back-transfer and rmse
    void setCheckQuality(bool x) { checkQuality = x; }
    // switch on/off weighted averaging/smoothing for cell data transfer
    // (default is off)
    void setContBool(bool x) { continuous = x; }
    // set the array names to name transferred data on target mesh
    void setNewArrayNames(const std::vector<std::string> &newnames) {
      newArrayNames = newnames;
    }
    // clear the new array names if set
    void unsetNewArrayNames() { newArrayNames.clear(); }

  protected:
    // number of points in mesh
    vtkIdType numPoints;
    // number of cells in mesh
    vtkIdType numCells;
    // set to true if sizeField attached to dataSet
    bool hasSizeField;
    // mesh points, topology and data
    vtkSmartPointer<vtkDataSet> dataSet;
    // name of mesh file
    std::string filename;
    // check transfer quality when on
    bool checkQuality;
    // switch on / off weighted averaging for cell data transfer
    // (default is off)
    bool continuous;
    // shape function order
    // (default is 1)
    int order;
    // new names to set for transferred data
    std::vector<std::string> newArrayNames;
    // --- for distributed data sets
    // --- (only populated for mesh resulting from call to meshBase::partition)
    // map between global and local node idx in partition
    std::map<vtkIdType, vtkIdType> globToPartNodeMap;
    // map between global and local cell idx in partition
    std::map<vtkIdType, vtkIdType> globToPartCellMap;
    // map between local and global node idx in partition
    std::map<vtkIdType, vtkIdType> partToGlobNodeMap;
    // map between local and global cell idx in partition
    std::map<vtkIdType, vtkIdType> partToGlobCellMap;
};

// --- auxiliary helpers
// compares two meshBase classes. used in testing
int diffMesh(meshBase *mesh1, meshBase *mesh2);

// sum comparison for vectors representing faces inserted into map
struct sortIntVec_compare {
  bool operator()(std::vector<int> lhs, std::vector<int> rhs) const;
};

// write patch map file for roc prep (trivial identity mapping)
void writePatchMap(const std::string &mapFile,
                   const std::map<int, int> &patchMap);
void writePatchMap(std::ofstream &outputStream,
                   const std::map<int, int> &patchMap);

#endif // meshBase_H
