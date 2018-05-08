#include <stdio.h>
#include <string.h>
#include <symmxGen.H>

// constructor
meshSymmx::meshSymmx(char* logFName, char* features, const char* licFName):
    haveLog(false),prog(NULL),model(NULL),mcase(NULL),mesh(NULL)
{
    // log
    if (strcmp(logFName,"NONE"))
    {
        haveLog = true;
        Sim_logOn(logFName);
    }
    // initialization
    SimLicense_start(features, licFName);
    MS_init();
    setProgress();
}

// destructor
meshSymmx::~meshSymmx()
{
    if (mesh)
        M_release(mesh);
    if (mcase)
        MS_deleteMeshCase(mcase);
    if (model)
        GM_release(model);
    SimLicense_stop();
}


void meshSymmx::setProgress()
{
    prog = Progress_new();
    Progress_setDefaultCallback(prog);    
}

void meshSymmx::createMeshFromModel(char* mFName)
{
  model = GM_load(mFName, 0, prog);
  mcase = MS_newMeshCase(model);
  MS_setMeshSize(mcase, GM_domain(model), 2, 0.5, 0);
  mesh = M_new(0, model);
  pSurfaceMesher surfaceMesher = SurfaceMesher_new(mcase, mesh);
  SurfaceMesher_execute(surfaceMesher, prog);
  SurfaceMesher_delete(surfaceMesher);
  pVolumeMesher volumeMesher = VolumeMesher_new(mcase, mesh);
  VolumeMesher_execute(volumeMesher, prog);
  VolumeMesher_delete(volumeMesher);
}

void meshSymmx::saveMesh(char* mFName)
{
    M_write(mesh, mFName, 0, prog);
}
