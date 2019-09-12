#ifndef CFMESHGEN_H
#define CFMESHGEN_H

#ifdef HAVE_CFMSH

#include "meshGen.H"

// OpenFOAM headers
#include <argList.H>
#include <fvOptions.H>

// VTK
#include <vtkUnstructuredGrid.h>

class cfmeshParams;

class cfmeshGen : public meshGen
{
  public:
    cfmeshGen();
    explicit cfmeshGen(cfmeshParams *params);
    ~cfmeshGen() override;

  // cfmesh mesh creation
  public:
    // create a 3D polygonal discretization from STL file
    int createMeshFromSTL(const char *fname) override;

    // reading openFoam mesh database
    void readFoamMesh();

  // internal management
  private:
    // pre-processor
    void initialize();

    // create dictionaries needed for openFoam
    void createControlDict();
    void createMshDict();
    void createfvSchemesDict();
    void createfvSolutionDict();

    // convert to built-in mesh database
    void genMshDB();
    void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                       int cellType, std::vector<int> &vrtIds);

    // surface edge extraction
    int surfaceFeatureEdgeDetect();

    // mesh quality improvement
    int improveMeshQuality();

  // class data and management flags
  private:
    bool defaults;
    std::string caseName;
    cfmeshParams *_params;

  // openfoam data structure
  private:
    Foam::argList *_args;
    Foam::Time *_runTime;
    Foam::fvMesh *_fmesh;
};

#endif

#endif
