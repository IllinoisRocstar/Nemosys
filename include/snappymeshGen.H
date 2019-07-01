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

class snappymeshGen : public meshGen
{
  // public
  public:
  snappymeshGen(); // Default constructor
  snappymeshGen(snappymeshParams* params); // Constructors w/ Parameters
  ~snappymeshGen(); // destructor

  public:
  // Creates mesh from input STL file
  int createMeshFromSTL(const char* fname);
  // Reads Foam mesh from polyMesh
  void readSnappyFoamMesh();
  // Writes mesh into polyMesh folder
  void writeMesh(const Foam::string& msg,
                const Foam::meshRefinement& meshRefiner,
                const Foam::meshRefinement::debugType debugLevel, 
                const Foam::meshRefinement::writeType writeLevel);

  private:
  // create dictionaries needed for openFoam
  void initialize();
  void createSnappyDict();
  void createControlDict();
  void createfvSchemesDict();
  void createfvSolutionDict();
  // Generates VTK mesh database for meshBase
  void genMshDB();
  void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet, 
  const int cellType, std::vector<int>& vrtIds);
  // Removes unneccesary patches from mesh
  void removeZeroSizedPatches(Foam::fvMesh& mesh);
  Foam::scalar getMergeDistance(const Foam::polyMesh& mesh,
                              const Foam::scalar mergeTol);
  // Extracts mesh surface from volumetric mesh
  void extractSurface(const Foam::polyMesh& mesh, const Foam::Time& runTime, 
  const Foam::labelHashSet& includePatches, const Foam::fileName& outFileName);
  // 
  Foam::label sizeCoeffToRefinement(const Foam::scalar level0Coeff,
                                  const Foam::scalar sizeCoeff);
  // Facilitates surface refinement over mesh
  Foam::autoPtr<Foam::refinementSurfaces> createRefinementSurfaces(
                            const Foam::searchableSurfaces& allGeometry,
                            const Foam::dictionary& surfacesDict,
                            const Foam::dictionary& shapeControlDict,
                            const Foam::label gapLevelIncrement,
                            const Foam::scalar level0Coeff);

  private:
  bool defaults;  // Default scenario
  std::string caseName; // Case name for internal use
  snappymeshParams* _params;  // SnappyHexMesh Parameters

  // OpenFOAM Data Structure
  private:
  Foam::argList* _args;
  Foam::Time* _runTime;
  Foam::fvMesh* _fmesh;
};

#endif

#endif