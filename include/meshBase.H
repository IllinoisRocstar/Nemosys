#ifndef meshBase_H
#define meshBase_H

// VTK
#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLHyperOctreeReader.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkHyperOctree.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtksys/SystemTools.hxx>
#include <vtkDoubleArray.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLHyperOctreeWriter.h>
#include <vtkXMLCompositeDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkCellCenters.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkGenericCell.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkCellLocator.h>
#include <vtkCellTreeLocator.h>
#include <vtkModifiedBSPTree.h>
#include <vtkMeshQuality.h>
#include <vtkQuadratureSchemeDefinition.h>
#include <vtkQuadratureSchemeDictionaryGenerator.h>
#include <vtkQuadraturePointsGenerator.h>
#include <vtkQuadraturePointInterpolator.h>
#include <vtkInformationQuadratureSchemeDefinitionVectorKey.h>
#include <baseInterp.H>
#include <spheres.H>
#include <vtkMath.h>

//MAd Lib
#include <MAdLib.h>
#include <ModelInterface.h>
#include <MeshDataBaseInterface.h>

// gmsh
#include <GmshEntities.h>
#include <Generator.h>

// netgen
namespace nglib
{
  #include<nglib.h>
}

// Nemosys
#include <AuxiliaryFunctions.H>

// stl
#include <functional>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <map>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <utility>

// TODO: 1) Implment copy and assignment constructors for all classes/subs
//       2) Migrate implementations to use smart pointer types 
// Virtual functions are implemented in vtkMesh.C

class meshingParams;

class meshBase 
{

  // constructors and destructors
  public:

    meshBase():numPoints(0),numCells(0),hasSizeField(0),checkQuality(0), order(1)
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
    // this method calls the other factor methods based on extension
    static meshBase* Create(std::string fname); 
    static std::shared_ptr<meshBase> CreateShared(std::string fname)
    {
      std::shared_ptr<meshBase> mesh;
      mesh.reset(meshBase::Create(fname));
      return mesh;
    }
    static std::unique_ptr<meshBase> CreateUnique(std::string fname)
    {
      return std::unique_ptr<meshBase>(meshBase::Create(fname));
    }
    // construct vtkMesh from gmsh msh file
    static meshBase* exportGmshToVtk(std::string fname);
    // construct vtkMesh from netgen vol file
    static meshBase* exportVolToVtk(std::string fname);
    // generate volume mesh from stl
    static meshBase* generateMesh(std::string fname, std::string meshEngine,
                                  meshingParams* params); 
    
  // access
  public:
    // get point with id
    virtual std::vector<double> getPoint(int id){}
    // get cell with id : returns point indices and respective coordinates
    virtual std::map<int, std::vector<double>> getCell(int id){}
    // get vector of coords of cell with id
    virtual std::vector<std::vector<double>> getCellVec(int id){}
    // get this meshes' dataSet
    vtkSmartPointer<vtkDataSet> getDataSet(){ return dataSet; }
    // register data to dataSet's point data
    virtual void setPointDataArray(const char* name, 
                                   const std::vector<std::vector<double>>& data){}
    // register data to dataSet's cell data
    virtual void setCellDataArray(const char* name,
                                  const std::vector<std::vector<double>>& data){}
    virtual void setCellDataArray(const char* name,
                                  const std::vector<double>& data){}
    // delete array with id from dataSet's point data
    virtual void unsetPointDataArray(int arrayID){}
    virtual void unsetPointDataArray(const char* name){}
    // delete array with id from dataSet's cell data
    virtual void unsetCellDataArray(int arrayID){} 
    virtual void unsetCellDataArray(const char* name){}
    // delete array with id from dataSet's field data
    virtual void unsetFieldDataArray(const char* name){}
    // get diameter of circumsphere of each cell
    virtual std::vector<double> getCellLengths(){}
    // get center of a cell
    virtual std::vector<double> getCellCenter(int cellID){}
    // build locators for efficient search operations
    vtkSmartPointer<vtkCellLocator> buildLocator();

    virtual void getIntegrationPointsAtCell(int cellID){}

  // data transfer
  public:
    // transfer point data with given id from this mesh to target
    // method can be "Finite Element", "Mortar Element", "RBF" etc.
    // only "Finite Element" has been implemented
    int transfer(meshBase* target, std::string method, 
                 const std::vector<int>& arrayIDs);
    int transfer(meshBase* target, std::string method, 
                 const std::vector<std::string>& arrayNames);
    // transfer all point data from this mesh to target
    int transfer(meshBase* target, std::string method);
 
  // size field generation
  public:
    // generate size field based on method ("gradient", "value", "error estimator")
    // and given a point data array. dev_mlt and maxIsmin are used to determine
    // which cells to consider for refinement
    void generateSizeField(std::string method, int arrayID, double dev_mlt, bool maxIsmin);
    // set whether this mesh has a sizefield
    void setSFBool(bool q) { hasSizeField = q; }
    // get this meshes' sizefield  
    bool getSFBool() { return hasSizeField; }
    // check for named array in vtk and return its integer id 
    int IsArrayName(std::string name);
		// set mesh order
		void setOrder(int _order) { order = _order;}
		int getOrder() { return order; } 
 
  // adaptive mesh refinement
  public:
    // perform sizefield-based h-refinement. edge_scale is for uniform refinement and 
    // is ignored in calls where method is "gradient", "value" etc. instead of "uniform"
    void refineMesh(std::string method, int arrayID, 
                    double dev_mult, bool maxIsmin, double edge_scale, 
										std::string ofname, bool transferData); 
    void refineMesh(std::string method, std::string arrayName, 
                    double dev_mult, bool maxIsmin, double edge_scale, 
										std::string ofname, bool transferData); 
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
    int getNumberOfPoints() {return numPoints;}
    // return the number of cells
    int getNumberOfCells() { return numCells;}
    void checkMesh(std::string ofname);


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
    std::string getFileName() { return filename; }
    void setCheckQuality(bool x) { checkQuality = x; }

  protected:
    int numPoints; // number of points in mesh
    int numCells;  // number of cells in mesh
    bool hasSizeField; // set to 1 if sizeField attached to dataSet
    vtkSmartPointer<vtkDataSet> dataSet; // mesh points, topology and data
    std::string filename; // name of mesh file
    bool checkQuality; // check transfer quality when on
		int order; // shape function order
};

// compares two meshBases- used in testing
int diffMesh(meshBase* mesh1, meshBase* mesh2);


#endif
