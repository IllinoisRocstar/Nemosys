// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>
#include <glob.h>

// Nemosys headers
#include <meshBase.H>
#include <meshPartitioner.H>
#include <meshStitcher.H>
#include <cgnsWriter.H>
#include <vtkAppendFilter.h>



std::vector<std::string> glob(const std::string& pattern)
{
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));
  int ret = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  if (ret == GLOB_NOSPACE || ret == GLOB_ABORTED)
  {
    globfree(&glob_result); 
    std::cerr << "glob() failed with return value " << ret << std::endl;
    exit(1);
  }
  else if (ret == GLOB_NOMATCH)
  {
    return std::vector<std::string>();
  }

  std::vector<std::string> fnames;
  for (int i = 0; i < glob_result.gl_pathc; ++i)
  {
    fnames.push_back(std::string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result); 
  return fnames;
}

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return glob(names.str());
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
  
  int surf = std::stoi(argv[4]);
  meshStitcher* stitcher = new meshStitcher(cgFnames, (bool) surf);
  delete stitcher; stitcher = nullptr;
  std::cout << "Application ended successfully!\n";
  return 0;
}

