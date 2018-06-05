#ifndef MESHSTITCHER_H
#define MESHSTITCHER_H

//std
#include <vector>
#include <string>

class cgnsAnalyzer;
class rocstarCgns;
class meshBase;

/* This class stitches together cgns grids and solution data. It exports the 
   stitched mesh with data to VTK data structures. The 1st component of 
   the member partitions, partitions[0], contains the stitched cgns object
   after instantiation of meshStitcher (when stitching vol/fluid files). 
    cgObj is the stitched mesh when stitching surf/ifluid files */
class meshStitcher
{
  public:
    // stitch fluid/vol files 
    meshStitcher(const std::vector<std::string>& cgFileNames);
    // stitch ifluid/surf files
    meshStitcher(int nCg, const std::string& baseCgName);
    meshStitcher(int begCg, int nCg, const std::string& baseCgName);
    ~meshStitcher();
    
    cgnsAnalyzer* getStitchedCGNS();
    meshBase* getStitchedMB();
 
  private:
    const std::vector<std::string> cgFileNames; // names of cgns files
    std::vector<cgnsAnalyzer*> partitions;      // vector of loaded cgns files
    meshBase* stitchedMesh;                     // stitched mesh in Nemosys format
    rocstarCgns* cgObj;
    // initialize series of cgns objects
    void initCgObj(int begCg, int nCg, const std::string& baseCgName);
    // disabling
    meshStitcher(const meshStitcher& that) = delete;
    meshStitcher& operator=(const meshStitcher& that) = delete;
    
};

#endif
