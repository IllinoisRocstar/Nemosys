#ifndef SYMMXGEN_H
#define SYMMXGEN_H

// std
#include <string>
// simmetrix
#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimModel.h"
#include "SimLicense.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include "SimDiscreteErrorCodes.h"
// Nemosys + vtk
#include <vtkUnstructuredGrid.h>
#include <meshGen.H>


// TODO: implement vtuTosymmx in-mem convertor that takes a vtkDataSet
//         as input 

class symmxParams;
    
class symmxGen : public meshGen
{

  public:
    // initialize params with default values
    symmxGen();
    symmxGen(symmxParams* params); 
    ~symmxGen();

  // symmetrix mesh creation
  public:
    // create mesh from symmetrix model file - (not currently used in nemo mesh gen)
    void createMeshFromModel(const char* mdlFName);
    // create model from stl file. relevant features (geomsim_discrete) must be passed
    // to object constructor to use this function
    int createModelFromSTL(const char* stlFName);
    int createSurfaceMeshFromSTL(const char* stlFName);
    int createVolumeMeshFromSTL(const char* stlFName);
    // base class mesh gen function
    int createMeshFromSTL(const char* fname); 
    // convert from symmetrix to vtu. if converting a volume mesh, you can choose to 
    // write only volume cells by setting writeSurfAndVol to false with the set method.
    void convertToVTU(); 
    void saveMesh(const std::string& mshFName);

  // access
  public:
    void setWriteSurfAndVol(bool b);

  // helpers
  private:
    void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                       const int numIds,
                       const int cellType,
                       pPList regionVerts);
    void addVtkVolCells(vtkSmartPointer<vtkUnstructuredGrid> dataSet);

  // status
  private:
    void setProgress();
    static void messageHandler(int type, const char* msg);

  private:
    symmxParams* params; // meshing params
    bool haveLog;
    bool writeSurfAndVol;
    pProgress prog;
    pGModel model;
    pDiscreteModel dModel;
    pACase mcase;
    pMesh mesh;
};

#endif
