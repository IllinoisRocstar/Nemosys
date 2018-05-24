#ifndef MESHSTITCHER_H
#define MESHSTITCHER_H

//std
#include <vector>
#include <string>

class cgnsAnalyzer;
class meshBase;

/* This class stitches together cgns grids and solution data. It exports the 
   stitched mesh with data to VTK data structures. The 1st component of 
   the member partitions, partitions[0], contains the stitched cgns object
   after instantiation of meshStitcher */
class meshStitcher
{
  public:
    
    meshStitcher(const std::vector<std::string>& cgFileNames);
    ~meshStitcher();
    
    cgnsAnalyzer* getStitchedCGNS();
    meshBase* getStitchedMB();
 
  private:
    const std::vector<std::string> cgFileNames; // names of cgns files
    std::vector<cgnsAnalyzer*> partitions;      // vector of loaded cgns files
    meshBase* stitchedMesh;                     // stitched mesh in Nemosys format
    // disabling
    meshStitcher(const meshStitcher& that) = delete;
    meshStitcher& operator=(const meshStitcher& that) = delete;
    
};

#endif
