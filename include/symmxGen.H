#ifndef SYMMXGEN_H
#define SYMMXGEN_H
#include <string>
#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimModel.h"
#include "SimLicense.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include "SimDiscreteErrorCodes.h"
// vtk
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

#include <vtkIdList.h>
#include <vtkCellTypes.h>

// method to write vtk grids
template<class TWriter> void writeVTFile(const std::string& fname, 
                                         vtkSmartPointer<vtkDataSet> dataSet)
{
  vtkSmartPointer<TWriter> Writer = vtkSmartPointer<TWriter>::New();
  Writer->SetFileName(&fname[0u]);
  Writer->SetInputData(dataSet);
  Writer->Write();
}

class meshSymmx
{

  public:
    /* possible feature options are :
      - geomsim_core
      - meshsim_surface
      - meshsim_volume
      - geomsim_discrete
      - meshsim_adv
      - meshsim_adapt */
    meshSymmx(const char* logFName="NONE", 
              const char* features="geomsim_core,meshsim_surface,meshsim_volume", 
              const char* licFName="simmodsuite.lic");
    ~meshSymmx();

  // symmetrix mesh creation
  public:
    // create mesh from symmetrix model file
    void createMeshFromModel(const char* mdlFName);
    // create model from stl file. relevant features (geomsim_discrete) must be passed
    // to object constructor to use this function
    int createModelFromSTL(const char* stlFName);
    int createMeshFromSTL(const char* stlFName);
    // convert from symmetrix to vtu 
    void convertToVTU();
    
    void saveMesh(const std::string& mshFName);


  // helpers
  private:
    void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                       const int numIds,
                       const int cellType,
                       pPList regionVerts);

  // status
  private:
    void setProgress();
    static void messageHandler(int type, const char* msg);

  private:
    bool haveLog;
    pProgress prog;
    pGModel model;
    pDiscreteModel dModel;
    pACase mcase;
    pMesh mesh;
    vtkSmartPointer<vtkUnstructuredGrid> dataSet;
};
#endif
