#ifndef MESHGENDRIVER_H
#define MESHGENDRIVER_H

#include <NemDriver.H>

class MeshGenDriver : public NemDriver
{
  public:
    
    MeshGenDriver(const std::string& ifname, const std::string& meshEngine, 
                  meshingParams* params, const std::string& ofname);
    ~MeshGenDriver();
    static MeshGenDriver* readJSON(const json& inputjson);
    static MeshGenDriver* readJSON(const std::string& ifname, const std::string& ofname, 
                                   const json& inputjson);
    meshBase* getNewMesh();
  private:
    meshBase* mesh; 
    meshingParams* params;
};

#endif
