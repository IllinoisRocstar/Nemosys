#ifndef SYMMXGEN_H
#define SYMMXGEN_H
#include <string>
#include "MeshSim.h"
#include "SimModel.h"
#include "SimLicense.h"
#include "SimError.h"
#include "SimErrorCodes.h"


class meshSymmx
{

  public:
    meshSymmx(char* logFName="NONE", char* features="geomsim_core,meshsim_surface,meshsim_volume", 
            const char* licFName="simmodsuite.lic");
    ~meshSymmx();

  // symmetrix mesh creation
  public:
    // create mesh from symmetrix model file
    void createMeshFromModel(char* mdlFName);
    void saveMesh(char* mshFName);

  private:
    void setProgress();

  private:
    bool haveLog;
    pProgress prog;
    pGModel model;
    pACase mcase;
    pMesh mesh;

};
#endif
