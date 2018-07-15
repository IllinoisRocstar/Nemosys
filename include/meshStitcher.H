#ifndef MESHSTITCHER_H
#define MESHSTITCHER_H

//std
#include <vector>
#include <string>
#include <memory>

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
    // stitch fluid/ifluid files depending on surf flag
    // surface stitching uses rocstarCgns, while vol uses cgnsAnalyzer 
    meshStitcher(const std::vector<std::string>& cgFileNames, bool surf);
    ~meshStitcher(){};
    
    std::shared_ptr<cgnsAnalyzer> getStitchedCGNS();
    std::shared_ptr<meshBase> getStitchedMB();
 
  private:
    // names of cgns files
    const std::vector<std::string> cgFileNames;             
    // vector of loaded cgns files
    std::vector<std::shared_ptr<cgnsAnalyzer>> partitions;
    // stitched mesh in Nemosys format  
    std::shared_ptr<meshBase> stitchedMesh;                     
    // stitched mesh in cgns fomrat
    std::shared_ptr<rocstarCgns> cgObj;
    // initialize series of cgns objects
    void initSurfCgObj();
    void initVolCgObj();
    // disabling
    meshStitcher(const meshStitcher& that) = delete;
    meshStitcher& operator=(const meshStitcher& that) = delete;
};

#endif
