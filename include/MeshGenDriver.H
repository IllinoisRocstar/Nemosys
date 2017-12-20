#ifndef MESHGENDRIVER_H
#define MESHGENDRIVER_H

#include <NemDriver.H>

class MeshGenDriver : public NemDriver
{
  public:
    
    MeshGenDriver(std::string ifname, std::string meshEngine, 
                  meshingParams* params, std::string ofname);
    ~MeshGenDriver();

    static MeshGenDriver* readJSON(json inputjson);

  private:
    meshBase* mesh; 
    meshingParams* params;
};

#endif
