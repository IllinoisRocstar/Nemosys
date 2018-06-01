#ifndef REMESH_DRIVER_H
#define REMESH_DRIVER_H

#include <NemDriver.H>

class meshStitcher;
class cgnsAnalyzer;

class RemeshDriver : public NemDriver
{
  public:
    RemeshDriver(const std::string& case_dir, const std::string& base_t,
                 const int fluidproc, const int ifluidniproc,
                 const int ifluidbproc, const int ifluidnbproc,
                 meshingParams* params);
    ~RemeshDriver();
    static RemeshDriver* readJSON(json inputjson);
  
  private:
    

    //std::vector<std::string> ifluidniNames; // ni cgns surf names
    //std::vector<std::string> ifluidnbNames; // nb cgns surf names
    //std::vector<std::string> ifluidbNames; // b cgns surf names
    std::vector<std::string> fluidNames; // vol cgns names
    std::string ifluidniName;
    std::string ifluidnbName;
    std::string ifluidbName;
    int fluidproc, ifluidniproc, ifluidbproc, ifluidnbproc;
    std::string base_t;
    std::string case_dir; 
    std::vector<meshStitcher*> stitchers; // stitches files in order listed above
    std::vector<cgnsAnalyzer*> cgObjs;  // stitched cgns in order listed above
    std::vector<meshBase*> mbObjs;   // stitched meshbase in order listed above
    meshingParams* params; // simmetrix meshing params 
    meshBase* stitchedSurf;
    meshBase* remeshedVol;   
    meshBase* remeshedSurf;
    void setCgFname(std::string& name, const std::string& prefix);
    void setCgFnames(std::vector<std::string>& names, const std::string& prefix, const int numproc);
    void symmxRemesh();
    void partitionMesh();
};

#endif
