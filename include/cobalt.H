#ifndef cobalt_H
#define cobalt_H

// stl
#include <memory>


class meshBase;

namespace COBALT
{

class cobalt
{

    public:
      cobalt(const std::shared_ptr<meshBase> fullMesh, const std::string inFnameVtk, const std::string outFnameCgr, const std::string outFnameCgi);

    public:
      void write();
      void writePatchMap(const std::string& mapFile, const std::map<int,int>& patchMap);
      void writePatchMap(std::ofstream& outputStream, const std::map<int,int>& patchMap);

    private:
      std::shared_ptr<meshBase> volMeshBase;
      std::shared_ptr<meshBase> surfMeshBase;
      std::string inFnameVtk;
      std::string outFnameCgr;
      std::string outFnameCgi;
};

}
#endif
