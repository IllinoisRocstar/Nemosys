#ifndef REMESH_DRIVER_H
#define REMESH_DRIVER_H

#include <NemDriver.H>

class meshStitcher;
class cgnsAnalyzer;
class MeshGenDriver;

class RemeshDriver : public NemDriver
{
  public:
    RemeshDriver(const std::string& case_dir, const std::string& base_t,
                 const int fluidproc, const int ifluidniproc,
                 const int ifluidbproc, const int ifluidnbproc, const json& remeshjson);
    ~RemeshDriver();
    static RemeshDriver* readJSON(json inputjson);
  
  private:
    // vol cgns names
    std::vector<std::string> fluidNames; 
    std::string ifluidniName;
    std::string ifluidnbName;
    std::string ifluidbName;
    int fluidproc, ifluidniproc, ifluidbproc, ifluidnbproc;
    std::string base_t;
    std::string case_dir; 
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
    void stitchFluidCGNS(const std::string& prefix, const int numproc);
    void stitchIFluidCGNS(std::string& name, const std::string& prefix, const int numproc);
    void setCgFname(std::string& name, const std::string& prefix);
    void setCgFnames(std::vector<std::string>& names, const std::string& prefix, const int numproc);
    void partitionMesh();
    // creates remeshed volume and surface
    void remesh(const json& remeshjson);
    // stitches the stitched cgns surfaces into one meshBase
    void stitchSurfaces();
};

#endif
