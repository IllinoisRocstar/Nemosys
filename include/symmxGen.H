#ifndef SYMMXGEN_H
#define SYMMXGEN_H
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
#include <meshGen.H>
//// vtk
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>


/* TODO: needs extension. currently, most parameters are hardcoded, except for
         the log file name, licence features and license file name/location 
    
      possible feature options are :
      - geomsim_core
      - meshsim_surface
      - meshsim_volume
      - geomsim_discrete
      - meshsim_adv
      - meshsim_adapt */

class SymmxParams : public meshingParams
{
  
  public:
    SymmxParams(){};

    ~SymmxParams(){};
  
  public:
    const char* logFName;
    const char* features;
    const char* licFName;
};


class symmxGen : public meshGen
{

  public:
    // initialize params with default values
    symmxGen();
    symmxGen(SymmxParams* params); 
    ~symmxGen();

  // symmetrix mesh creation
  public:
    // create mesh from symmetrix model file
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
    //vtkSmartPointer<vtkDataSet> getDataSet();
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
    bool haveLog;
    bool writeSurfAndVol;
    pProgress prog;
    pGModel model;
    pDiscreteModel dModel;
    pACase mcase;
    pMesh mesh;
//    vtkSmartPointer<vtkDataSet> dataSet;
};
#endif
