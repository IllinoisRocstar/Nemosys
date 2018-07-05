#ifndef meshBase_H
#define meshBase_H

// VTK
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

//MAd Lib
#include <MAdLib.h>

// Nemosys
#include <pntMesh.H>

// stl
#include <vector>
#include <iosfwd>
#include <fstream>
#include <map>
#include <string>
#include <memory>


// NOTE: virtual methods are usually implemented in vtkMesh.C

class meshingParams;

class meshBase 
{

  // constructors and destructors
  public:

    meshBase()
      : dataSet(0),numPoints(0),numCells(0),
        hasSizeField(0),checkQuality(0), continuous(1),order(1)
    {
      std::cout << "meshBase constructed" << std::endl;
    }

    virtual ~meshBase()
    {
      std::cout << "meshBase destroyed" << std::endl;
    }
  
  // meshBase factory methods
  public:
    // construct vtkMesh from filename
    // this method calls the other factory methods based on extension
    // caller must delete object after use
    static meshBase* Create(std::string fname); 
    // construct from existing vtkDataSet and assign newname as filename
    // caller must delete object after use
    static meshBase* Create(vtkSmartPointer<vtkDataSet> other, std::string newname);
    // create shared ptr from fname
    // memory is managed by shared pointer, so do not call delete after use
    static std::shared_ptr<meshBase> CreateShared(std::string fname);
    // create unique ptr from fname
    // memory is managed by unique pointer, so do not call delete after use
    static std::unique_ptr<meshBase> CreateUnique(std::string fname);
    /* create from coordinates and connectivities.
       NOTE: use of this is only valid when mesh has one cell type. 
             cellType one of the vtkCellType enums.
             Currently, only VTK_TETRA and VTK_TRIANGLE are supported
       Caller must delete object after use  */
    static meshBase* Create(const std::vector<double>& xCrds,
                            const std::vector<double>& yCrds,
                            const std::vector<double>& zCrds,
                            const std::vector<int>& elmConn, const int cellType,
                            std::string newname);

    // construct vtkMesh from gmsh msh file (called in Create methods)
    static meshBase* exportGmshToVtk(std::string fname);
    // construct vtkMesh from netgen vol file (called in Create methods)
    static meshBase* exportVolToVtk(std::string fname);
    // construct vtkMesh from netgen vol file (called in Create methods)
    static meshBase* exportPntToVtk(std::string fname);
  
    // generate volume mesh from stl
    // caller must delete object after use
    static meshBase* generateMesh(std::string fname, std::string meshEngine,
                                  meshingParams* params); 
    // stitch together several meshBases
    // caller must delete object after use
    static meshBase* stitchMB(const std::vector<meshBase*>& mbObjs);
    // mesh partitioning (with METIS)
    // caller must delete objects after use
    static std::vector<meshBase*> partition(const meshBase* mbObj, const int numPartitions);
      
 
  // access
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
    virtual void inspectEdges(const std::string& ofname){}
    // get this meshes' dataSet
    vtkSmartPointer<vtkDataSet> getDataSet() const { return dataSet; }
    // extract the surface mesh
    virtual vtkSmartPointer<vtkDataSet> extractSurface(){}

    // register data to dataSet's point data
    virtual void setPointDataArray(const char* name, 
                                   const std::vector<std::vector<double>>& data){}
    virtual void setPointDataArray(const char* name, const std::vector<double>& data){}
    // register data to dataSet's cell data
    virtual void setCellDataArray(const char* name,
                                  const std::vector<std::vector<double>>& data){}
    virtual void setCellDataArray(const char* name,
                                  const std::vector<double>& data){}
    // get scalar point or cell data array. assumes data is not allocated prior to calling
    virtual void getPointDataArray(std::string name,
                                   std::vector<double>& data) {}
    virtual void getCellDataArray(std::string name,
                                  std::vector<double>& data) {}
    // delete array with id from dataSet's point data
    virtual void unsetPointDataArray(int arrayID){}
    virtual void unsetPointDataArray(const char* name){}
    // delete array with id from dataSet's cell data
    virtual void unsetCellDataArray(int arrayID){} 
    virtual void unsetCellDataArray(const char* name){}
    // delete array with id from dataSet's field data
    virtual void unsetFieldDataArray(const char* name){}
    // get diameter of circumsphere of each cell
    virtual std::vector<double> getCellLengths() const {}
    // get center of a cell
    virtual std::vector<double> getCellCenter(int cellID) const {}
    // build locators for efficient search operations
    vtkSmartPointer<vtkCellLocator> buildLocator();
    // get cell type as an integer
    // assumes all elements are the same type
    virtual int getCellType() const = 0;
    virtual std::vector<int> getConnectivities() const {}

  // data transfer
  public:
    // transfer point data (pointOrCell=0) or cell data (pointOrCell=1)
    //  with given id from this mesh to target
    // method can be "Consistent Interpolation", "Mortar Element", "RBF" etc.
    // only "Consistent Interpolation" has been implemented
    int transfer(meshBase* target, std::string method, 
                 const std::vector<int>& arrayIDs, const bool pointOrCell = 1);
    int transfer(meshBase* target, std::string method, 
                 const std::vector<std::string>& arrayNames, const bool pointOrCell = 1);
    // transfer all point and cell data from this mesh to target
    int transfer(meshBase* target, std::string method);

  // integration
  public:
    // integrate arrays in arrayIDs over the mesh. 
    // integrated data is available per cell after operation
    // total integral for each datum is returned
    std::vector<std::vector<double>> 
    integrateOverMesh(const std::vector<int>& arrayIDs);

 
  // size field generation
  public:
    // generate size field based on method ("gradient", "value", "error estimator")
    // and given a point data array. dev_mlt and maxIsmin are used to determine
    // which cells to consider for refinement
    void generateSizeField(std::string method, int arrayID, double dev_mlt, 
                           bool maxIsmin, double sizeFactor=1.);
    // set whether this mesh has a sizefield
    void setSFBool(bool q) { hasSizeField = q; }
    // get this meshes' sizefield  
    bool getSFBool() const { return hasSizeField; }
    // check for named array in vtk and return its integer id (0-point data, 1-cell data)
    int IsArrayName(std::string name, const bool pointOrCell = 0);
    // set mesh order
    void setOrder(int _order) { order = _order;}
    int getOrder() const { return order; } 
 
  // adaptive mesh refinement
  public:
    // perform sizefield-based h-refinement. edge_scale is for uniform refinement and 
    // is ignored in calls where method is "gradient", "value" etc. instead of "uniform"
    void refineMesh(std::string method, int arrayID, 
                    double dev_mult, bool maxIsmin, double edge_scale, 
                    std::string ofname, bool transferData, double sizeFactor=1.); 
    void refineMesh(std::string method, std::string arrayName, 
                    double dev_mult, bool maxIsmin, double edge_scale, 
                    std::string ofname, bool transferData, double sizeFactor=1.); 
    // added for uniform refinement by driver
    void refineMesh(std::string method, double edge_scale, std::string ofname, bool transferData); 
    void refineMesh(std::string method, int arrayID, int order, 
                    std::string ofname, bool transferData);
    void refineMesh(std::string method, std::string arrayName, int order, 
                    std::string ofname, bool transferData);
  // diagnostics
  public:
    // generate a report of the mesh
    virtual void report(){}; 
    // return the number of points
    int getNumberOfPoints() const {return numPoints;}
    // return the number of cells
    int getNumberOfCells() const { return numCells;}
    void checkMesh(std::string ofname);

  // for distributed data sets
  public:
    std::map<int,int> getGlobToPartNodeMap() { return globToPartNodeMap; }
    std::map<int,int> getGlobToPartCellMap() { return globToPartCellMap; }
    std::map<int,int> getPartToGlobNodeMap() { return partToGlobNodeMap; }
    std::map<int,int> getPartToGlobCellMap() { return partToGlobCellMap; }
  
  // write and conversion
  public:
    // write the mesh to file named after the private var 'filename'
    virtual void write(){}
    // write the mesh to file named fname
    virtual void write(std::string fname){}
    // convert to gmsh format without data 
    void writeMSH(std::ofstream& outputStream);
    void writeMSH(std::string fname);
    // convert to gmsh format with specified point or cell data
    void writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID); 
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID);
    // convert to gmsh format with specified point or cell data for
    // only volume elements (USE ONLY FOR MADLIB STUFF)
    void writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID, 
                  bool onlyVol); // added for overloading, doesn't do anything 
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID,
                  bool onlyVol);

    void setFileName(std::string fname) { filename = fname; }
    std::string getFileName() const { return filename; }
    void setCheckQuality(bool x) { checkQuality = x; }
    void setContBool(bool x) { continuous = x;}
    void setNewArrayNames(const std::vector<std::string>& newnames) { newArrayNames = newnames; }
    void unsetNewArrayNames() { newArrayNames.clear(); }

  protected:
    int numPoints; // number of points in mesh
    int numCells;  // number of cells in mesh
    bool hasSizeField; // set to 1 if sizeField attached to dataSet
    vtkSmartPointer<vtkDataSet> dataSet; // mesh points, topology and data
    std::string filename; // name of mesh file
    bool checkQuality; // check transfer quality when on
    bool continuous; // switch on / off weighted averaging for cell data transfer
    int order; // shape function order
    std::vector<std::string> newArrayNames;

    // for distributed data sets
    std::map<int,int> globToPartNodeMap; // map between global and local node idx in partition
    std::map<int,int> globToPartCellMap; // map between global and local cell idx in partition
    std::map<int,int> partToGlobNodeMap;
    std::map<int,int> partToGlobCellMap;
};

// compares two meshBases- used in testing
int diffMesh(meshBase* mesh1, meshBase* mesh2);

#endif
