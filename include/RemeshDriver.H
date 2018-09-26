#ifndef REMESH_DRIVER_H
#define REMESH_DRIVER_H

#include <NemDriver.H>
#include <memory>

class meshStitcher;
class cgnsAnalyzer;
class MeshGenDriver;

class RemeshDriver : public NemDriver
{
  public:
    RemeshDriver(const std::vector<std::string>& fluidNames,
                 const std::vector<std::string>& ifluidniNames,
                 const std::vector<std::string>& ifluidnbNames,
                 const std::vector<std::string>& ifluidbNames,
                 const std::vector<std::string>& _burnNames,
                 const std::vector<std::string>& _iBurnNames,
                 const json& remeshjson, int numPartitions,
                 const std::string& base_t, int writeIntermediateFiles, 
                 double searchTolerance, const std::string& caseName);
    ~RemeshDriver();
    static RemeshDriver* readJSON(json inputjson);
 
  private:
    // cgns file names
    std::vector<std::string> fluidNames;
    std::vector<std::string> ifluidniNames;
    std::vector<std::string> ifluidnbNames;
    std::vector<std::string> ifluidbNames; 
    std::vector<std::string> burnNames;
    std::vector<std::string> iBurnNames;
    // mesh generator which uses remeshing parameters from json
    std::unique_ptr<MeshGenDriver> mshgendrvr;  
    // stitches files in order listed above
    std::vector<std::unique_ptr<meshStitcher>> stitchers;
    // stitched cgns in order listed above  
    std::vector<std::shared_ptr<cgnsAnalyzer>> cgObjs;
    // stitched meshbase in order listed above    
    std::vector<std::shared_ptr<meshBase>> mbObjs;  
    // all surface cgns files stitched together 
    std::shared_ptr<meshBase> stitchedSurf;
    // surface burn files stitched together
    std::shared_ptr<meshBase> stitchedBurnSurf;

    // stitched volume cgns remeshed
    std::shared_ptr<meshBase> remeshedVol;  
    // surface extracted from remeshed stitched volume 
    std::shared_ptr<meshBase> remeshedSurf;
    // number of partitions
    int numPartitions;

  //helpers
  private:
    // stitch surf/vol cgns files, using rocstarCgns if surf=1
    void stitchCGNS(const std::vector<std::string>& fnames, bool surf);
    // creates remeshed volume and surface
    void remesh(const json& remeshjson);
    // stitches the stitched cgns surfaces into one meshBase
    void stitchSurfaces();
};

// get cg file names using rocstar convention in regex glob
inline std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                            const std::string& prefix,
                                            const std::string& base_t);

#endif
