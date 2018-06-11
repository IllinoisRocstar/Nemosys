#ifndef REMESH_DRIVER_H
#define REMESH_DRIVER_H

#include <NemDriver.H>

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
                 const json& remeshjson);
    ~RemeshDriver();
    static RemeshDriver* readJSON(json inputjson);
    void writeCobalt(const std::string& mapFile, const std::string& ofname);
 
  private:
    // cgns file names
    std::vector<std::string> fluidNames;
    std::vector<std::string> ifluidniNames;
    std::vector<std::string> ifluidnbNames;
    std::vector<std::string> ifluidbNames; 
    // mesh generator
    MeshGenDriver* mshgendrvr;  
    // stitches files in order listed above
    std::vector<meshStitcher*> stitchers;
    // stitched cgns in order listed above  
    std::vector<cgnsAnalyzer*> cgObjs;
    // stitched meshbase in order listed above    
    std::vector<meshBase*> mbObjs;  
    // all surface cgns files stitched together 
    meshBase* stitchedSurf;
    // stitched volume cgns remeshed
    meshBase* remeshedVol;  
    // surface extracted from remeshed stitched volume 
    meshBase* remeshedSurf;
  
  //helpers
  private:
    // stitch surf/vol cgns files, using rocstarCgns if surf=1
    void stitchCGNS(const std::vector<std::string>& fnames, bool surf);
    // write patch map file for roc prep (trivial identity mapping)
    void writePatchMap(const std::string& mapFile, const std::map<int,int>& patchMap); 
    void writePatchMap(std::ofstream& outputStream, const std::map<int,int>& patchMap); 
    void writeCobalt(const std::string& mapFile, std::ofstream& outputStream);
    // partitions mesh (NOT USED)
    void partitionMesh();
    // creates remeshed volume and surface
    void remesh(const json& remeshjson);
    // stitches the stitched cgns surfaces into one meshBase
    void stitchSurfaces();

};

// get cg file names using rocstar convention in glob
inline std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                            const std::string& prefix,
                                            const std::string& base_t);

// sum comparison for vectors representing faces inserted into map
struct sortIntVec_compare
{
  bool operator() (std::vector<int> lhs, std::vector<int> rhs) const;
};

#endif
