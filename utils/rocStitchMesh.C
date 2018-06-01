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

void setCgFnames(std::vector<std::string>& names, const std::string& prefix,
                 const std::string& base_t, const int numproc)
{
  names.resize(numproc);
  for (int i = 0; i < numproc; ++i)
  {
    std::stringstream basename;
    basename << prefix << "_" << base_t << "_";
    if (i < 10)
      basename << "000" << i << ".cgns";
    else if (i >= 10 && i < 100)
      basename << "00" << i << ".cgns";
    else if (i >= 100 && i < 1000)
      basename << 0 << i << ".cgns";
    names[i] = basename.str();
  } 
}

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  int nInCgFile = 0;
  
  std::vector<std::string> cgFileName;
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " nCgFile CgFileName0 surf?" << std::endl
              << "Eg) rocStitchMesh 4 fluid_04.124000_0000.cgns 0" << std::endl; 
    return 0;
  }
  std::string::size_type sz;   // alias of size_t
  nInCgFile = std::stoi(argv[1],&sz);
  int surf = std::stoi(argv[3]);
  if (surf)
  {
    meshStitcher* stitcher = new meshStitcher(nInCgFile, argv[2]);
    delete stitcher; stitcher = 0;
  }
  else
  {
    std::string base_t(argv[2]);
    std::size_t pos = base_t.find_first_of("_");
    base_t = base_t.substr(pos+1,9);
    std::vector<std::string> fluidNames;
    setCgFnames(fluidNames, "fluid", base_t, nInCgFile);
    meshStitcher* stitcher = new meshStitcher(fluidNames);
    delete stitcher; stitcher = 0;
  }     
 

  std::cout << "Application ended successfully!\n";
  return 0;
}



