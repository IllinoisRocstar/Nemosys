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

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  int nInCgFile = 0;
  
  std::vector<std::string> cgFileName;
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " nCGNSFileIn inCgFileName1 inCgFileName2 ..." << std::endl;
    return 0;
  }
  std::string::size_type sz;   // alias of size_t
  nInCgFile = std::stoi(argv[1],&sz);
  for (int iCg=0; iCg < nInCgFile; iCg++)
     cgFileName.push_back(argv[2+iCg]);

  meshStitcher* stitcher = new meshStitcher(cgFileName);
  meshBase* trgVTK = stitcher->getStitchedMB();
 
  if (stitcher){ delete stitcher; trgVTK = 0;}
  
  std::cout << "Application ended successfully!\n";

  return 0;
}
