#ifndef patran_H
#define patran_H

#include "nemosys_export.h"
#include "meshBase.H"

namespace PATRAN {

class NEMOSYS_EXPORT patran
{
  public:
    patran(std::shared_ptr<meshBase> fullMesh,
           const std::string &inFnameVtk,
           const std::string &outFnameNeu,
           const std::map<int, int> &faceTypeMap,
           const std::map<int, int> &nodeTypeMap,
           const std::map<int, bool> &nodeStructuralMap,
           const std::map<int, bool> &nodeMeshMotionMap,
           const std::map<int, bool> &nodeThermalMap,
           const std::vector<int> &nppItr);

  public:
    void write25(std::ofstream &outputStream) const;
    void write26(std::ofstream &outputStream) const;
    void write1(std::ofstream &outputStream) const;
    void write2(std::ofstream &outputStream) const;
    void write6(std::ofstream &outputStream);
    void write8(std::ofstream &outputStream);
    void write99(std::ofstream &outputStream) const;
    bool comparePatch(int i, int j);

  private:
    std::shared_ptr<meshBase> fullMesh;
    std::string inFnameVtk;
    std::string outFnameNeu;
    std::shared_ptr<meshBase> volMeshBase;
    std::shared_ptr<meshBase> surfMeshBase;
    std::map<int, int> faceTypeMap;
    std::map<int, int> nodeTypeMap;
    std::map<int, bool> nodeStructuralMap;
    std::map<int, bool> nodeMeshMotionMap;
    std::map<int, bool> nodeThermalMap;
    std::vector<int> nppVec;
    std::map<int, std::string> face2nodes;
    std::map<nemId_t, std::vector<int>> boundaryNodeId2PatchNo;

};

} // namespace PATRAN

#endif
