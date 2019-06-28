#ifndef cobalt_H
#define cobalt_H

#include "nemosys_export.h"
#include "meshBase.H"

namespace COBALT {

class NEMOSYS_EXPORT cobalt
{
  public:
    cobalt(std::shared_ptr<meshBase> fullMesh,
           const std::string &inFnameVtk,
           const std::string &outFnameCgr,
           const std::string &outFnameCgi);

  public:
    void write() const;
    void writePatchMap(const std::string &mapFile,
                       const std::map<nemId_t, nemId_t> &patchMap) const;
    void writePatchMap(std::ofstream &outputStream,
                       const std::map<nemId_t, nemId_t> &patchMap) const;

  private:
    std::shared_ptr<meshBase> volMeshBase;
    std::shared_ptr<meshBase> surfMeshBase;
    std::string inFnameVtk;
    std::string outFnameCgr;
    std::string outFnameCgi;
};

} // namespace COBALT

#endif
