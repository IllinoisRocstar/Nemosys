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
#ifdef HAVE_GLOB_H
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
#else
  return nemAux::glob(case_dir, prefix, base_t);
#endif
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

  int surf = std::stoi(argv[4]);
  bool is_surf = (bool) surf;

  // remove surface file names in case still in the volume file names
  if (!is_surf) {
    std::vector<std::string> newCgFnames;
    for (auto it = cgFnames.begin(); it != cgFnames.end(); it++) {
      if (it->find("i" + std::string(argv[2])) == std::string::npos)
        newCgFnames.push_back(*it);
    }
    cgFnames = newCgFnames;
  }

  nemAux::printVec(cgFnames);
  meshStitcher *stitcher = new meshStitcher(cgFnames, is_surf);
  delete stitcher; stitcher = nullptr;
  std::cout << "Application ended successfully!\n";
  return 0;
}

