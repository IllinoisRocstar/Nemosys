/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
// standard headers
#include <cstring>
#include <sstream>
#include <iostream>
#include <memory>

// Nemosys headers
#include <Mesh/meshBase.H>
#include <MeshPartitioning/meshPartitioner.H>
#include <MeshPartitioning/meshStitcher.H>
#include <IO/cgnsWriter.H>
#include <vtkAppendFilter.h>
#ifdef HAVE_GLOB_H
#  include <glob.h>
#  include <cstring>
#else
#  include <boost/filesystem.hpp>
#endif

namespace {

// check if _exp is in the _name string
bool expCont(const std::string &_name, const std::string &_exp) {
  return (_name.find(_exp) != std::string::npos);
}

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
    if (expCont(_timestamp, "*")) overide = true;

    std::string fext(it->path().filename().string());
    if (!expCont(fext, _prefix)) continue;
    {
      if (fext[0] != _prefix[0]) continue;

      if (!overide)
        if (!expCont(fext, _timestamp)) continue;

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

  for (const auto &cgFname : cgFnames) {
    std::cout << cgFname << " ";
  }
  std::cout << '\n';
  meshStitcher *stitcher = new meshStitcher(cgFnames, is_surf);
  delete stitcher; stitcher = nullptr;
  std::cout << "Application ended successfully!\n";
  return 0;
}

