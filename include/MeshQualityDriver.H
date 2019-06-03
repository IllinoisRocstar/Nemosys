#ifndef MESHQUALITYDRIVER_H
#define MESHQUALITYDRIVER_H

#include <NemDriver.H>

class MeshQualityDriver : public NemDriver
{
  public:

    MeshQualityDriver() : mesh(NULL) {};
    MeshQualityDriver(std::string _mesh, std::string ofname); 
    ~MeshQualityDriver();

    static MeshQualityDriver* readJSON(json inputjson);

  private:
    meshBase* mesh; 
};

#endif
