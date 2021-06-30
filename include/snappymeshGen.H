#ifndef SNAPPYMESHGEN_H
#define SNAPPYMESHGEN_H

#ifdef HAVE_CFMSH 

#include "meshGen.H"
#include "argList.H"
#include "fvOptions.H"
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>
#include <vtkCellArray.h> 

// snappyHexMesh Headers
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *// 

class snappymeshParams;


/** 
    @brief snappymeshGen facilitates full-hexahedral/hex-dominent meshing of 
          complex geometries with surface, patch, or region refinement 
          capabilities in addition with adding layers to geometry. A 
          background mesh is needed in order to perform meshing operation.
          Inputs can be given in form of surface files (.stl, .obj, etc.).
          Output mesh is written in OpenFOAM polyMesh format.
          snappymeshGen is based on snappyHexMesh utility of OpenFOAM with
          following features
          - Castallated Mesh
          - Snap
          - Add Layers 
**/
class snappymeshGen : public meshGen
{
  // --- constructors and destructor
  public:
  /** @brief snappymeshGen standard constructor
  **/
  snappymeshGen(); // Default constructor

  /** @brief snappymeshGen alternate constructor
      @param params snappymeshParams object
  **/
  snappymeshGen(snappymeshParams* params); // Constructors w/ Parameters

  /** @brief snappymeshGen standard destructor
  **/
  ~snappymeshGen(); // destructor

  // --- access
  public:
  /** @brief Creates mesh from input STL file
  **/
  int createMeshFromSTL(const char* fname) override;

  /** @brief Reads Foam mesh from polyMesh
  **/
  void readSnappyFoamMesh();
  
  /** @brief Writes mesh into polyMesh folder
  **/
  void writeMesh(const Foam::string& msg,
                const Foam::meshRefinement& meshRefiner,
                const Foam::meshRefinement::debugType debugLevel, 
                const Foam::meshRefinement::writeType writeLevel);

  private:  
  /** @brief Initializes OpenFOAM environment
  **/
  void initialize();
  
  /** @brief Creates snappyHexMeshDict for meshing operation
  **/
  void createSnappyDict();
  
  /** @brief Creates dictionary for OpenFOAM runtime environment
  **/
  void createControlDict();
  
  /** @brief Creates dictionary for OpenFOAM runtime environment
  **/
  void createfvSchemesDict();
  
  /** @brief Creates dictionary for OpenFOAM runtime environment
  **/
  void createfvSolutionDict();
  
  /** @brief Generates meshBase database for Foam mesh after mesh writing 
            operation
  **/
  void genMshDB();
  
  /** @brief Creates and adds VTK cell to incoming database using celltype and
            cell bounding point Ids
      @param dataset VTK unstructured dataset
      @param cellType VTK cell type
      @param vrtIds Cell bounding point Ids
  **/
  void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet, 
  const int cellType, std::vector<int>& vrtIds);
  
  /** @brief Removes unneccesary patches from mesh
  **/
  void removeZeroSizedPatches(Foam::fvMesh& mesh);
  
  /** @brief
  **/
  Foam::scalar getMergeDistance(const Foam::polyMesh& mesh,
                              const Foam::scalar mergeTol);
  
  /** @brief Extracts mesh surface from volumetric mesh
  **/
  void extractSurface(const Foam::polyMesh& mesh, const Foam::Time& runTime, 
  const Foam::labelHashSet& includePatches, const Foam::fileName& outFileName);
  
  /** @brief Calculated refinement ration based on base level coefficient and 
            size coefficient
  **/
  // 
  Foam::label sizeCoeffToRefinement(const Foam::scalar level0Coeff,
                                  const Foam::scalar sizeCoeff);
  
  /** @brief Facilitates surface refinement over mesh
  **/
  Foam::autoPtr<Foam::refinementSurfaces> createRefinementSurfaces(
                            const Foam::searchableSurfaces& allGeometry,
                            const Foam::dictionary& surfacesDict,
                            const Foam::dictionary& shapeControlDict,
                            const Foam::label gapLevelIncrement,
                            const Foam::scalar level0Coeff);

  private:
  /** @brief If enabled, generated mesh using default parameters
  **/
  bool defaults;  // Default scenario

  /** @brief snappymeshParams object Parameters
  **/
  snappymeshParams* _params;  // SnappyHexMesh Parameters

  // -- OpenFOAM Data Structure
  private:
  Foam::argList* _args;
  Foam::Time* _runTime;
  Foam::fvMesh* _fmesh;
};

#endif

#endif