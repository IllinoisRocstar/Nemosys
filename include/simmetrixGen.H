#ifndef SIMMETRIXGEN_H
#define SIMMETRIXGEN_H

// Nemosys headers
#include "nemosys_export.h"
#include "meshGen.H"

// simmetrix
#include <MeshSim.h>

// VTK
#include <vtkUnstructuredGrid.h>
#include <ModelTypes.h>

// TODO: implement vtuToSimmetrix in-mem converter that takes a vtkDataSet
//       as input

class simmetrixParams;

class NEMOSYS_EXPORT simmetrixGen : public meshGen
{
  public:
    // initialize params with default values
    simmetrixGen();
    explicit simmetrixGen(simmetrixParams *params);
    ~simmetrixGen() override;

  // simmetrix mesh creation
  public:
    // create mesh from simmetrix model file
    // (not currently used in Nemosys mesh gen)
    void createMeshFromModel(const char *mdlFName);
    // create model from stl file. relevant features (geomsim_discrete) must be
    // passed to object constructor to use this function
    int createModelFromSTL(const char *stlFName);
    int createSurfaceMeshFromSTL(const char *stlFName);
    int createVolumeMeshFromSTL(const char *stlFName);
    // base class mesh gen function
    int createMeshFromSTL(const char *fname) override;
    // convert from simmetrix to vtu. if converting a volume mesh, you can
    // choose to write only volume cells by setting writeSurfAndVol to false
    // with the set method.
    void convertToVTU();
    void saveMesh(const std::string &mshFName);

  // access
  public:
    void setWriteSurfAndVol(bool b);

  // helpers
  private:
    void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                       int numIds,
                       int cellType,
                       pPList regionVerts);
    void addVtkVolCells(vtkSmartPointer<vtkUnstructuredGrid> dataSet);

  // status
  private:
    void setProgress();
    static void messageHandler(int type, const char *msg);

  private:
    simmetrixParams *params; // meshing params
    bool haveLog;
    bool writeSurfAndVol;
    pProgress prog;
    pGModel model;
    pDiscreteModel dModel;
    pACase mcase;
    pMesh mesh;
};

#endif
