%module pyNemosysSimmetrix
%import "pyNemosys.i"
%{
#include "RocRestartDriver.H"
#include "simmetrixGen.H"
#include "simmetrixParams.H"
#include "cgnsAnalyzer.H"
#include "meshStitcher.H"
#include "meshPartitioner.H"
%}


class simmetrixParams : public meshingParams {
  public:
    simmetrixParams();
    ~simmetrixParams();

    std::string logFName;
    std::string features;
    std::string licFName;
    double meshSize;
    double anisoMeshCurv;
    double minCurvSize;
    double glbSizeGradRate;
    double surfMshImprovGradRate;
    double surfMshImprovMinSize;
};


class simmetrixGen : public meshGen {
  public:
    simmetrixGen();
    simmetrixGen(simmetrixParams *params);
    ~simmetrixGen();

    void createMeshFromModel(const char *mdlFName);
    int createModelFromSTL(const char *stlFName);
    int createSurfaceMeshFromSTL(const char *stlFName);
    int createVolumeMeshFromSTL(const char *stlFName);
    int createMeshFromSTL(const char *fname);
    void convertToVTU();
    void saveMesh(const std::string &mshFName);
    void setWriteSurfAndVol(bool b);
};


class meshStitcher {
  public:
    meshStitcher(const std::vector<std::string> &cgFileNames, bool surf);
    ~meshStitcher();

    std::shared_ptr<cgnsAnalyzer> getStitchedCGNS();
    std::shared_ptr<meshBase> getStitchedMB();
};


class meshPartitioner {
  public:
    meshPartitioner(int nNde, int nElm, std::vector<int> &elemConn,
                    MeshType_t meshType);
    meshPartitioner(cgnsAnalyzer *inCg);
    meshPartitioner(meshBase *inMB);
    // destructor
    ~meshPartitioner();

    // mesh information
    int partition(int nPartition);
    int partition();
    std::vector<double> getCrds(int iPart, std::vector<double> crds);
    std::vector<int> getConns(int iPart);
    int getNNdePart(int iPart);
    int getNElmPart(int iPart);
};
