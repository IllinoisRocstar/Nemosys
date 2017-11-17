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
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkCellCenters.h>
#include <vtkCell.h>
#include <baseInterp.H>
#include <spheres.H>

//MAd Lib
#include <MAdLib.h>
#include <ModelInterface.h>
#include <MeshDataBaseInterface.h>

// gmsh
#include <GmshEntities.h>
#include <Generator.h>

// other
#include<netgenInterface.H>
#include <AuxiliaryFunctions.H>

#include <StlToVtk.H>

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
    // construct vtkMesh surface from STL
    static meshBase* exportStlToVtk(std::string fname); 
  
  // access
  public:
    // get point with id
    virtual std::vector<double> getPoint(int id){};
    // get cell with id : returns point indices and respective coordinates
    virtual std::map<int, std::vector<double>> getCell(int id){};
   
  // diagnostics
  public:
    virtual void report(char* fname){}; 
    int getNumberOfPoints() {return numPoints;}
    int getNumberOfCells() { return numCells;}
    virtual void write(std::string fname){};



  protected:
    int numPoints;
    int numCells;
    vtkDataSet* dataSet; // must manually delete
};

class meshUser
{
  public:
    meshUser(std::string _fname): mesh(),fname(_fname)
    {
      mesh = meshBase::Create(fname);
      std::cout << "user constructed" << std::endl;
    }
    
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
  
    // diagnostics
    public:
      // print coordinates of point with id
      void printPoint(int id);
      // print point ids and coordinates of cell with id
      void printCell(int id);

  private:
    meshBase* mesh;
    std::string fname;

};

// Auxiliary Functions
template<typename T>
void printVec(const std::vector<T>& v);
//std::string trim_fname(std::string fname, std::string ext);


#endif
