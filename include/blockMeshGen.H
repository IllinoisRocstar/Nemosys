#ifndef BLOCKMESHGEN_H
#define BLOCKMESHGEN_H

#ifdef HAVE_CFMSH

// Nemosys Headers
#include "meshGen.H"
#include "foamMesh.H"
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>
#include <vtkCellArray.h>

// Foam headers
#include "argList.H"
#include "fvOptions.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

/*
  blockMeshGen <-- meshGen <-- meshBase
  This class incorporates mesh generating method of blockMesh utility.
  Currently, we are supporting automatic generation and full-hex meshing
  of Box, Sphere, Cylinder, and Tapered Cone geometries.
*/

class blockMeshParams;  // Friend Class

class blockMeshGen : public meshGen
{
  public:
  blockMeshGen();   // Default constructor
  blockMeshGen(blockMeshParams* params);  // Constructor with parameters
  ~blockMeshGen();  // Class destructor


  public:
  // Generates mesh and returns VTK database
  int createMeshFromSTL(const char* fname);
  // Reads OpenFOAM mesh from polyMesh
  void readFoamMesh();


  private:
  // Initializes OpenFOAM args and runtime
  void initialize();
  // Creates blockMeshDict from user arguments
  void createBlockMshDict();
  // Creates dummy controlDict
  void createControlDict();
  // Creates dummy fvSchemesDict
  void createfvSchemesDict();
  // Creates dummy fvSolutionDict
  void createfvSolutionDict();
  // Generated VTK database from Foam mesh. Used in readFoamMesh function
  void genMshDB();
  // Generated VTK mesh statistics
  void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                    const int cellType, std::vector<int>& vrtIds);

  private:
  bool defaults;  // Default parameters boolean
  std::string caseName; // Meshing case name for internal use
  blockMeshParams* _params; // Definition of parameters pointer

  // OpenFOAM Data Structure
  private:
  Foam::argList* _args;
  Foam::Time* _runTime;
  Foam::fvMesh* _fmesh;
  
};


#endif

#endif
