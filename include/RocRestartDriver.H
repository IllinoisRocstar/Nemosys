#ifndef ROCRESTARTDRIVER_H
#define ROCRESTARTDRIVER_H

#include <NemDriver.H>

class meshStitcher;
class cgnsAnalyzer;
 
class RocRestartDriver : public NemDriver
{

  public:
    RocRestartDriver(const std::vector<std::string>& fluidNamesRm,
                     const std::vector<std::string>& ifluidniNamesRm,
                     const std::vector<std::string>& ifluidnbNamesRm,
                     const std::vector<std::string>& ifluidbNamesRm,
                     const std::vector<std::string>& fluidNamesLts,
                     const std::vector<std::string>& ifluidniNamesLts,  
                     const std::vector<std::string>& ifluidnbNamesLts,   
                     const std::vector<std::string>& ifluidbNamesLts); 
 
    ~RocRestartDriver(); 
    static RocRestartDriver* readJSON(json inputjson);

  private:
    // remeshed cgns file names
    std::vector<std::string> fluidNamesRm;
    std::vector<std::string> ifluidniNamesRm;
    std::vector<std::string> ifluidnbNamesRm;
    std::vector<std::string> ifluidbNamesRm;
    // last ts cgns file names
    std::vector<std::string> fluidNamesLts;
    std::vector<std::string> ifluidniNamesLts;
    std::vector<std::string> ifluidnbNamesLts;
    std::vector<std::string> ifluidbNamesLts;
    // stitches files in order listed above
    std::vector<meshStitcher*> stitchers;
    // stitched last ts meshbase in order listed above    
    std::vector<meshBase*> mbObjs;  
    // all last ts surface cgns files stitched together (mbObjs[1:-1])
    meshBase* stitchedSurf;
    // remeshed, 0 ts loaded cgns files
    std::vector<cgnsAnalyzer*> fluidRmCg;
    std::vector<cgnsAnalyzer*> ifluidNiRmCg;
    std::vector<cgnsAnalyzer*> ifluidNbRmCg;
    std::vector<cgnsAnalyzer*> ifluidBRmCg;  
    // remeshed, 0 ts loaded cgns files converted to vtu
    std::vector<meshBase*> fluidRmMb;
    std::vector<meshBase*> ifluidNiRmMb;
    std::vector<meshBase*> ifluidNbRmMb;
    std::vector<meshBase*> ifluidBRmMb;  


  //helpers
  private:
    // populates stitchers
    void stitchCGNS(const std::vector<std::string>& fnames, bool surf);
    // populates cgns and meshbase vectors
    void loadPartCgMb();
    std::pair< std::vector<cgnsAnalyzer*>, std::vector<meshBase*> >
      loadCGNS(const std::vector<std::string>& fnames, bool surf);
    // stitch ni, b and nb surfaces
    void stitchSurfaces();
    // transfer stitched data to new meshbase partitions
    void transferStitchedToPartMb();

};

// get cg file names using rocstar convention in glob
inline std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                            const std::string& prefix,
                                            const std::string& base_t);
#endif
