// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include <meshBase.H>
#include <meshPartitioner.H>
#include <meshStitcher.H>
#include <cgnsWriter.H>
#include <vtkAppendFilter.h>
#include <AuxiliaryFunctions.H>

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
}

int main(int argc, char* argv[])
{
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " case_dir prefix base_t surf=0|1" << std::endl
              << "Eg) rocStitchMesh ~/my_case fluid 04.124000 0" << std::endl; 
    return 0;
  }
  
  std::vector<std::string> cgFnames(getCgFNames(argv[1], argv[2], argv[3]));
  printVec(cgFnames);

  cgFnames.pop_back();
  cgFnames.pop_back();
  cgFnames.pop_back();
  
  int surf = std::stoi(argv[4]);
  meshStitcher* stitcher = new meshStitcher(cgFnames, (bool) surf);
  delete stitcher; stitcher = nullptr;
  std::cout << "Application ended successfully!\n";
  return 0;
}

