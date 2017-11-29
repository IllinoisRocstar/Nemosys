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
#include <vtkTetra.h>
#include <baseInterp.H>
#include <spheres.H>

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
//#include <meshGen.H>

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

/*TODO: This class is responsible for reading in a *mesh* file and exporting it 
        to a vtkMesh object. In this way, we have access to the methods provided
          through the vtkDataSet- the core member of a vtkMesh object.
  
  TODO: MESHGEN: Given a cad file, surface or volume mesh generation will be delegated to 
        a separate mesh gen class

  TODO: TRANSFER: Given two vtkMesh objects of same geometry with point or cell data, 
        transfer data from one to other. This will be plugged in after the refinement
        and will be based on MadLibs or implementing using methods from VTK API
  
  TODO: SFGEN: This class will provide the interface for generating size fields based
        on data-> derived classes like gradient, value, error estimator etc
        output is a vtkMesh object with size field data
  
*/

class meshBase
{

  // constructors and destructors
  public:

    meshBase():numPoints(0),numCells(0),dataSet(NULL)
    {
      std::cout << "meshBase constructed" << std::endl;
    }
    virtual ~meshBase()
    {
      if(dataSet) dataSet->Delete();
      std::cout << "meshBase destroyed" << std::endl;
    }
  
  // factory methods (kind of like 'virtual constructors')
  // these functions return objects derived from meshBase
  public:
    // construct vtkMesh from filename
    // this method calls the other factor methods based on extension
    static meshBase* Create(std::string fname); 
    // construct vtkMesh from gmsh msh file
    static meshBase* exportGmshToVtk(std::string fname);
    // construct vtkMesh from netgen vol file
    static meshBase* exportVolToVtk(std::string fname);
    // generate volume mesh from stl
    static meshBase* generateMesh(std::string fname, std::string meshEngine);    

    
  // access
  public:
    // get point with id
    virtual std::vector<double> getPoint(int id){};
    // get cell with id : returns point indices and respective coordinates
    virtual std::map<int, std::vector<double>> getCell(int id){};
    vtkDataSet* getDataSet(){ return dataSet; }
	  // returns cells sharing point (topological search)
	  virtual std::vector<int> getCellsWithPoint(int pnt){};
    virtual void setPointDataArray(const char* name, 
              std::vector<std::vector<double>>& data){};
    virtual void unsetPointDataArray(int arrayID){};                       
  // data transfer
  public:
    // transfer data from this mesh to target
    int transfer(meshBase* target, std::string method);
 
  // diagnostics
  public:
    
    virtual void report(char* fname){}; 
    int getNumberOfPoints() {return numPoints;}
    int getNumberOfCells() { return numCells;}

  //write
  public:
    virtual void write(std::string fname, std::string extension){};

  protected:
    int numPoints;
    int numCells;
    vtkDataSet* dataSet; // must manually delete
};

class meshUser
{
  public:
    
    meshUser():mesh(){};

    meshUser(std::string _fname): mesh(),fname(_fname)
    {
      mesh = meshBase::Create(fname);
      std::string ext = find_ext(fname);
      // setting write extension for proper vtk file type output
      write_ext = ext; 
      if (ext == ".vtk" || ext == ".vol" || ext == ".msh") // default is vtu
      {
        std::string new_ext(".vtu");
        write_ext = new_ext;
      }  
      
      std::cout << "user constructed" << std::endl;
    }

    int generateMesh(std::string filename, std::string meshEngine); //meshEngine=netgen||gmsh
    
    ~meshUser()
    {
      if (mesh)
      {
        delete mesh;
        mesh = 0;
      }
      std::cout << "user destroyed" << std::endl;
    }

    // information
    public:
      // get number of points in mesh
      int getNumberOfPoints(); 
      // get number of cells in mesh
      int getNumberOfCells();
      // generate a report of the mesh
      void report();      

    // access
    public:
      // get point with id
      std::vector<double> getPoint(int id);
      // get cell with id : returns point indices and respective coordinates
      std::map<int,std::vector<double>> getCell(int id); 
      // get mesh object
      meshBase* getMesh(){ return mesh; };

  
    // data transfer
    public:
      int transfer(meshUser* target, std::string method);


    // diagnostics
    public:
      // print coordinates of point with id
      void printPoint(int id);
      // print point ids and coordinates of cell with id
      void printCell(int id);

    // write
    public:
      void write(); 
      void write(std::string fname);      
  private:
    meshBase* mesh;
    std::string fname;
    std::string write_ext;
};


#endif
