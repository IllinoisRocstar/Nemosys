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
                 const std::string& base_t, const int numproc, const int beg)
{
  names.resize(numproc);
  int j = 0;
  for (int i = beg; i < numproc+beg; ++i)
  {
    std::stringstream basename;
    basename << prefix << "_" << base_t << "_";
    if (i < 10)
      basename << "000" << i << ".cgns";
    else if (i >= 10 && i < 100)
      basename << "00" << i << ".cgns";
    else if (i >= 100 && i < 1000)
      basename << 0 << i << ".cgns";
    names[j] = basename.str();
    j+=1;
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
              << " nCgFile CgFileName0 part_num surf?" << std::endl
              << "Eg) rocStitchMesh 4 fluid_04.124000_0001.cgns 1 0" << std::endl; 
    return 0;
  }
  std::string::size_type sz;
  nInCgFile = std::stoi(argv[1],&sz);
  int surf = std::stoi(argv[4]);
  int part_num = std::stoi(argv[3]);
  if (surf)
  {
    meshStitcher* stitcher = new meshStitcher(part_num, nInCgFile, argv[2]);
    delete stitcher; stitcher = 0;
  }
  else
  {
    std::string base_t(argv[2]);
    std::size_t pos = base_t.find_first_of("_");
    std::string prefix = base_t.substr(0,pos);
    base_t = base_t.substr(pos+1,9);
    std::vector<std::string> fluidNames;
    setCgFnames(fluidNames, prefix, base_t, nInCgFile, part_num);
    printVec(fluidNames);
    meshStitcher* stitcher = new meshStitcher(fluidNames);
    delete stitcher; stitcher = 0;
  }     
 

  std::cout << "Application ended successfully!\n";
  return 0;
}



