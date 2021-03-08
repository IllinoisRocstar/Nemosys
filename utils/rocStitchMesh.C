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
#ifdef HAVE_GLOB_H
#  include <glob.h>
#  include <cstring>
#else
#  include <boost/filesystem.hpp>
#endif

namespace {
// search for pattern and return vector of matches
#ifdef HAVE_GLOB_H
std::vector<std::string> glob(const std::string &pattern) {
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));
  int ret = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);

  if (ret == 0) {
    std::vector<std::string> fnames;
    for (std::size_t i = 0; i < glob_result.gl_pathc; ++i)
      fnames.emplace_back(glob_result.gl_pathv[i]);
    globfree(&glob_result);
    return fnames;
  } else if (ret == GLOB_NOMATCH) {
    globfree(&glob_result);
    return std::vector<std::string>();
  } else {  // if (ret == GLOB_NOSPACE || ret == GLOB_ABORTED) {
    std::cerr << "glob() failed with return value " << ret << std::endl;
    globfree(&glob_result);
    exit(1);
  }
}
#else
std::vector<std::string> glob(const std::string &_src_dir,
                              const std::string &_prefix,
                              const std::string &_timestamp) {
  std::vector<std::string> fnames;
  std::cout << _src_dir << " " << _prefix << " " << _timestamp << std::endl;

  // recursive copy all content that contain given name
  // iterate and delete matches
  for (boost::filesystem::directory_iterator end_dir_it, it(_src_dir);
       it != end_dir_it; ++it) {
    bool overide = false;
    if (nemAux::expCont(_timestamp, "*")) overide = true;

    std::string fext(it->path().filename().string());
    if (!nemAux::expCont(fext, _prefix)) continue;
    {
      if (fext[0] != _prefix[0]) continue;

      if (!overide)
        if (!nemAux::expCont(fext, _timestamp)) continue;

      std::cout << "Identified " << it->path() << std::endl;
      fnames.push_back(it->path().filename().string());
      // boost::filesystem::copy_file(it->path(), des_path,
      // boost::filesystem::copy_option::overwrite_if_exists);
    }
  }
  return fnames;
}
#endif
}  // namespace

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
#ifdef HAVE_GLOB_H
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return glob(names.str());
#else
  return glob(case_dir, prefix, base_t);
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

