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

class blockMeshParams;  // Friend Class


/**
  @brief blockMeshGen <-- meshGen <-- meshBase
        This class incorporates mesh generating method of blockMesh utility.
        Currently, we are supporting automatic full-hex mesh generation
        Box, Sphere, Cylinder, and Tapered Cone geometries. Output mesh is
        written in OpenFOAM polyMesh format.
**/
class blockMeshGen : public meshGen
{
  // --- constructors and destructor
  public:
  /** @brief blockMeshGen standard constructor
  **/
  blockMeshGen();   // Default constructor

  /** @brief blockMeshGen alternate constructor. Uses user-defined parameters
            to perform requested meshing operation.
      @param params blockMeshParams object
  **/
  blockMeshGen(blockMeshParams* params);  // Constructor with parameters

  /** @brief blockMeshGen standard desctructor
  **/
  ~blockMeshGen();  // Class destructor

  // -- meshGen methods override
  public:
  /** @brief Generates mesh and returns VTK database
      @param fname Input surface file name
      @return meshBase mesh dataset
  **/
  int createMeshFromSTL(const char* fname);

  /** @brief Reads OpenFOAM mesh from polyMesh directory
  **/
  void readFoamMesh();


  // -- Internal
  private:
  /**
      @brief Initializes OpenFOAM args and runtime
  **/
  void initialize();

  /**
      @brief Creates blockMeshDict from user arguments
  **/
  void createBlockMshDict();

  /**
      @brief Creates dummy controlDict
  **/
  void createControlDict();

  /**
      @brief Creates dummy fvSchemesDict
  **/
  void createfvSchemesDict();

  /**
      @brief Creates dummy fvSolutionDict
  **/
  void createfvSolutionDict();

  /**
      @brief Generated VTK database from Foam mesh. Used in readFoamMesh function
  **/
  void genMshDB();

  /**
      @brief Generated VTK mesh statistics
  **/
  void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                    const int cellType, std::vector<int>& vrtIds);

  // --- Internal Use
  private:
  bool defaults;  // Default parameters boolean

  /** @brief Definition of parameters pointer
  **/
  blockMeshParams* _params; // Definition of parameters pointer

  // --- OpenFOAM Data Structure
  private:
  Foam::argList* _args{};
  Foam::Time* _runTime{};
  Foam::fvMesh* _fmesh{};

};


#endif

#endif
