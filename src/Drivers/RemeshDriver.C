#ifdef HAVE_SIMMETRIX
// Nemosys headers
#  include "AuxiliaryFunctions.H"
#  include "RemeshDriver.H"
#  include "TransferDriver.H"

#  include "MeshGenDriver.H"
#  include "RocPartCommGenDriver.H"
#  include "cgnsAnalyzer.H"
#  include "meshStitcher.H"

#  include <utility>

// file system
#  include <boost/filesystem.hpp>

namespace NEM {
namespace DRV {

void backUpDir(const std::string &_cur_dir, const std::string &_bak_dir);
void rmDirContent(const boost::filesystem::path &_dir, bool _skipDir);
bool expContains(const std::string &_name, const std::string &_exp);
void rmFileName(const boost::filesystem::path &_dir, const std::string &_fn);
void cpFileName(const boost::filesystem::path &_src_dir,
                const boost::filesystem::path &_des_dir, const std::string &_fn,
                const std::string &_starts);
bool createDir(const boost::filesystem::path &directory);

RemeshDriver::RemeshDriver(
    std::vector<std::string> _fluidNames,
    std::vector<std::string> _ifluidniNames,
    std::vector<std::string> _ifluidnbNames,
    std::vector<std::string> _ifluidbNames, std::vector<std::string> _burnNames,
    std::vector<std::string> _iBurnNames, const jsoncons::json &remeshjson,
    int _numPartitions, const std::string &base_t, int writeIntermediateFiles,
    double searchTolerance, const std::string &caseName,
    const std::map<std::string, std::vector<int>> &surfacePatchTypes,
    const std::string &_prefix_path)
    : fluidNames(std::move(_fluidNames)),
      ifluidniNames(std::move(_ifluidniNames)),
      ifluidnbNames(std::move(_ifluidnbNames)),
      ifluidbNames(std::move(_ifluidbNames)),
      burnNames(std::move(_burnNames)),
      iBurnNames(std::move(_iBurnNames)),
      numPartitions(_numPartitions) {
  // setting path
  prefixPath = (_prefix_path.length() != 0) ? _prefix_path + "/" : "./";

  // preparing input for the actual operator class
  // stitch fluid files
  this->stitchCGNS(fluidNames, false);
  // stitch burn files
  this->stitchCGNS(burnNames, false);
  // stitch iburn files
  this->stitchCGNS(iBurnNames, true);
  // stitch ifluid_ni files
  this->stitchCGNS(ifluidniNames, true);
  // stitch ifluid_b files
  this->stitchCGNS(ifluidbNames, true);
  // stitch ifluid_nb files
  this->stitchCGNS(ifluidnbNames, true);

  // get stitched meshes from stitcher vector
  for (int i = 0; i < this->stitchers.size(); ++i) {
    // cgObjs.push_back(stitchers[i]->getStitchedCGNS());
    this->mbObjs.push_back(stitchers[i]->getStitchedMB());
    this->mbObjs[i]->setContBool(false);
  }
  if (writeIntermediateFiles)
    this->mbObjs[0]->write(prefixPath + "stitchedVol.vtu");
  // creates remeshedVol and remeshedSurf
  this->remesh(remeshjson);
  // creates stitchedSurf
  if (this->mbObjs.size() > 1) {
    // CGNS files are already stitched together into meshBase obj
    // need to transfer this data onto the surface mesh

    // stitches b, ni and nb surfaces into one stitchedSurf
    // also stitches iburn files
    this->stitchSurfaces();
    // do not smooth for cell data transfer
    this->stitchedSurf->setContBool(false);
    this->stitchedBurnSurf->setContBool(false);
    // transfer patch number to remeshed surface
    std::vector<std::string> transferPaneData{"patchNo", "bcflag",
                                              "cnstr_type"};
    /*
    this->stitchedSurf->transfer(remeshedSurf.get(), "Consistent Interpolation",
                                 transferPaneData, true);
                                 */
    auto transfer = TransferDriver::CreateTransferObject(
        this->stitchedSurf.get(), remeshedSurf.get(),
        "Consistent Interpolation");
    transfer->transferCellData(transferPaneData);

    if (writeIntermediateFiles) this->stitchedSurf->write();
    if (writeIntermediateFiles) this->stitchedBurnSurf->write();
  }

  if (writeIntermediateFiles) this->remeshedSurf->write();
  bool withC2CTransSmooth = remeshjson.contains("C2CTransSmooth")
                                ? remeshjson["C2CTransSmooth"].as<bool>()
                                : false;

  // instantiating and executing actual operator class
  std::unique_ptr<RocPartCommGenDriver> rocprepdrvr =
      std::unique_ptr<RocPartCommGenDriver>(new RocPartCommGenDriver(
          this->remeshedVol, this->remeshedSurf, this->mbObjs[0],
          this->stitchedSurf, this->stitchedBurnSurf, numPartitions, base_t,
          writeIntermediateFiles, searchTolerance, caseName, surfacePatchTypes,
          withC2CTransSmooth, prefixPath));
  std::cout << "RemeshDriver created" << std::endl;
}

RemeshDriver::~RemeshDriver() {
  std::cout << "RemeshDriver destroyed" << std::endl;
}

void RemeshDriver::stitchCGNS(const std::vector<std::string> &fnames,
                              bool surf) {
  if (!fnames.empty()) {
    stitchers.push_back(
        std::unique_ptr<meshStitcher>(new meshStitcher(fnames, surf)));
  }
}

void RemeshDriver::remesh(const jsoncons::json &remeshjson) {
  std::cout << "Extracting surface mesh" << std::endl;
  std::unique_ptr<meshBase> surf = std::unique_ptr<meshBase>(meshBase::Create(
      mbObjs[0]->extractSurface(), prefixPath + "extractedSurface.vtp"));
  // remeshing with engine specified in input
  surf->write(prefixPath + "extractedSurface.stl");
  mshgendrvr = std::unique_ptr<MeshGenDriver>(
      MeshGenDriver::readJSON(prefixPath + "extractedSurface.stl",
                              prefixPath + "remeshedVol.vtu", remeshjson));
  // get the remeshed volume from the mesh generator
  remeshedVol = mshgendrvr->getNewMesh();
  // extract the remeshed volumes surface
  remeshedSurf = meshBase::CreateShared(remeshedVol->extractSurface(),
                                        prefixPath + "remeshedSurf.vtp");
}

void RemeshDriver::stitchSurfaces() {
  // stitch b, ni and nb surfaces
  std::vector<std::shared_ptr<meshBase>> surfs;
  surfs.insert(surfs.begin(), mbObjs.begin() + 3, mbObjs.end());
  stitchedSurf = meshBase::stitchMB(surfs);
  stitchedSurf->setFileName(prefixPath + "stitchedSurf.vtu");

  std::vector<std::shared_ptr<meshBase>> burnSurfs;
  burnSurfs.insert(burnSurfs.begin(), mbObjs.begin() + 2, mbObjs.begin() + 3);
  stitchedBurnSurf = meshBase::stitchMB(burnSurfs);
  stitchedBurnSurf->setFileName(prefixPath + "stitchedBurnSurf.vtu");
}

RemeshDriver *RemeshDriver::readJSON(const jsoncons::json &inputjson) {
  std::string case_name = inputjson["Rocstar Case Name"].as<std::string>();
  std::string case_dir = inputjson["RocFluMP Case Directory"].as<std::string>();
  std::string burn_dir =
      inputjson["RocBurnAPN Case Directory"].as<std::string>();
  std::string base_t = inputjson["Base Time Step"].as<std::string>();
  std::vector<std::string> fluNames(getCgFNames(case_dir, "fluid", base_t));
  std::vector<std::string> ifluniNames(
      getCgFNames(case_dir, "ifluid_ni", base_t));
  std::vector<std::string> iflunbNames(
      getCgFNames(case_dir, "ifluid_nb", base_t));
  std::vector<std::string> iflubNames(
      getCgFNames(case_dir, "ifluid_b", base_t));
  std::vector<std::string> burnNames(getCgFNames(burn_dir, "burn", base_t));
  std::vector<std::string> iBurnNames(getCgFNames(burn_dir, "iburn", base_t));
  if (fluNames.empty() || burnNames.empty()) {
    std::cerr
        << "Error finding Rocstar output files. Make sure folder locations are"
           " properly set in the input file and all folders are non-empty."
        << std::endl;
    throw;
  }
  jsoncons::json remeshjson = inputjson["Remeshing Options"];
  int numPartitions = inputjson["Number Of New Partitions"].as<int>();
  int writeIntermediateFiles = 0;
  double searchTolerance;
  if (inputjson.contains("Search Tolerance"))
    searchTolerance = inputjson["Search Tolerance"].as<double>();
  else
    searchTolerance = 1.e-16;

  if (inputjson.contains("Write Intermediate Files")) {
    writeIntermediateFiles = inputjson["Write Intermediate Files"].as<int>();
  }

  // Map surface patch numbers to surface types
  std::map<std::string, std::vector<int>> surfacePatchTypes;

  if (inputjson.contains("Burning Surface Patches")) {
    std::vector<int> bPatchMap =
        inputjson["Burning Surface Patches"].as<std::vector<int>>();
    // sort by patch ID
    std::sort(bPatchMap.begin(), bPatchMap.end());
    if (!bPatchMap.empty()) {
      surfacePatchTypes["Burning"] = bPatchMap;
    }
  } else {
    std::cerr << "Specify Burning Surface Patches in input file." << std::endl;
    throw;
  }
  if (inputjson.contains("Non-Burning Surface Patches")) {
    std::vector<int> nbPatchMap =
        inputjson["Non-Burning Surface Patches"].as<std::vector<int>>();
    // sort by patch ID
    std::sort(nbPatchMap.begin(), nbPatchMap.end());
    if (!nbPatchMap.empty()) {
      surfacePatchTypes["Non-Burning"] = nbPatchMap;
    }
  } else {
    std::cerr << "Specify Non-Burning Surface Patches in input file."
              << std::endl;
    throw;
  }
  if (inputjson.contains("Non-Interacting Surface Patches")) {
    std::vector<int> niPatchMap =
        inputjson["Non-Interacting Surface Patches"].as<std::vector<int>>();
    // sort by patch ID
    std::sort(niPatchMap.begin(), niPatchMap.end());
    if (!niPatchMap.empty()) {
      surfacePatchTypes["Non-Interacting"] = niPatchMap;
    }
  } else {
    std::cerr << "Specify Non-Interacting Surface Patches in input file."
              << std::endl;
    throw;
  }

  int writeToTemporaryDirectory = 0;
  std::string tmp_main_dir;
  if (inputjson.contains("Write To Temporary Directory")) {
    writeToTemporaryDirectory =
        inputjson["Write To Temporary Directory"].as<int>();
    if (writeToTemporaryDirectory) {
      if (inputjson.contains("Temporary Rocstar Case Directory")) {
        tmp_main_dir =
            inputjson["Temporary Rocstar Case Directory"].as<std::string>();
      } else {
        std::cerr << "Specify Temporary RocfluMP Case Directory and"
                     " Temporary RocburnAPN Case Directory in input file."
                  << std::endl;
        throw;
      }
    }
  }
  // set temporary dir
  const char dir_path[] = "/remesh_tmp";

  // preforming prep actions
  boost::filesystem::path cwd_dir;
  cwd_dir = boost::filesystem::current_path();
  boost::filesystem::path rmsh_dir(cwd_dir);
  rmsh_dir += boost::filesystem::path(dir_path);
  std::cout << "Current working directory " << cwd_dir << std::endl;

  // make a tmp directory or clean if already there
  if (boost::filesystem::is_directory(rmsh_dir)) {
    std::cout << "Removing existing " << rmsh_dir << std::endl;
    int nRem = boost::filesystem::remove_all(rmsh_dir);
    std::cout << "Removed " << nRem << " item(s)." << std::endl;
    if (boost::filesystem::create_directory(rmsh_dir))
      std::cout << "Created " << rmsh_dir << std::endl;
    else {
      std::cerr << "Problem creating directory " << rmsh_dir << std::endl;
      throw;
    }
  } else if (boost::filesystem::create_directory(rmsh_dir)) {
    std::cout << "Created " << rmsh_dir << std::endl;
  }

  // remeshing, repartitioning and solution transfer actions
  RemeshDriver *remeshdrvobj = new RemeshDriver(
      fluNames, ifluniNames, iflunbNames, iflubNames, burnNames, iBurnNames,
      remeshjson, numPartitions, base_t, writeIntermediateFiles,
      searchTolerance, case_name, surfacePatchTypes, rmsh_dir.string());

  if (writeToTemporaryDirectory) {
    // backing up all old cgns files for rocflu (TODO: expensive remove)
    boost::filesystem::remove_all(tmp_main_dir);
    backUpDir(case_dir + "/../..", tmp_main_dir);

    // set output file to temporary directory
    case_dir = tmp_main_dir + "/Rocflu/Rocout";
    burn_dir = tmp_main_dir + "/RocburnAPN/Rocout";
  }

  // remove rocout files from remesh time step
  rmFileName(case_dir, base_t);
  rmFileName(burn_dir, base_t);

  // backing up all old cgns files for rocflu (TODO: expensive remove)
  backUpDir(case_dir, rmsh_dir.string() + "/_flu_bak/");
  backUpDir(burn_dir, rmsh_dir.string() + "/_burn_bak/");

  // post action
  // remove old cgns files (TODO: copy somewhere)
  // remove old cgns and files
  rmFileName(case_dir, ".cgns");
  rmFileName(case_dir, ".txt");
  rmFileName(burn_dir, ".cgns");
  rmFileName(burn_dir, ".txt");

  // remove flu-modin files
  boost::filesystem::path flumodin_dir(case_dir + "/../Modin/");
  rmFileName(flumodin_dir, ".map");
  rmFileName(flumodin_dir, ".dim");
  rmFileName(flumodin_dir, ".com");
  rmFileName(flumodin_dir, ".cmp");
  rmFileName(flumodin_dir, ".rnm");

  // copy over old cgns and txt files
  cpFileName(rmsh_dir.string() + "/_flu_bak/", case_dir, "fluid", "fluid");
  cpFileName(rmsh_dir.string() + "/_flu_bak/", case_dir, "ifluid", "ifluid");
  cpFileName(rmsh_dir.string() + "/_burn_bak/", burn_dir, "burn", "burn");
  cpFileName(rmsh_dir.string() + "/_burn_bak/", burn_dir, "iburn", "iburn");

  // copy over newly generated cgns and txt files
  cpFileName(rmsh_dir, case_dir, "fluid", "fluid");
  cpFileName(rmsh_dir, case_dir, "ifluid", "ifluid");
  cpFileName(rmsh_dir, burn_dir, "burn", "burn");
  cpFileName(rmsh_dir, burn_dir, "iburn", "iburn");

  // copying newly generated rocflu modin files
  cpFileName(rmsh_dir, flumodin_dir, ".map", "*");
  cpFileName(rmsh_dir, flumodin_dir, ".dim", "*");
  cpFileName(rmsh_dir, flumodin_dir, ".com", "*");
  cpFileName(rmsh_dir, flumodin_dir, ".cmp", "*");

  return remeshdrvobj;
}

std::vector<std::string> getCgFNames(const std::string &case_dir,
                                     const std::string &prefix,
                                     const std::string &base_t) {
  std::stringstream names;
#  ifdef HAVE_GLOB_H
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
#  else
  return nemAux::glob(case_dir, prefix, base_t);
#  endif
}

// file transfer helpers
// recursive copy, will create destination folder (make sure it is removed
// first)
bool copyDir(const boost::filesystem::path &source,
             const boost::filesystem::path &destination) {
  namespace fs = boost::filesystem;
  try {
    // Check whether the function call is valid
    if (!fs::exists(source) || !fs::is_directory(source)) {
      std::cerr << "Source directory " << source.string()
                << " does not exist or is not a directory." << std::endl;
      return false;
    }
    if (fs::exists(destination)) {
      std::cerr << "Destination directory " << destination.string()
                << " already exists." << std::endl;
      return false;
    }
    // Create the destination directory
    if (!fs::create_directory(destination)) {
      std::cerr << "Unable to create destination directory"
                << destination.string() << std::endl;
      return false;
    }
  } catch (fs::filesystem_error const &e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
  // Iterate through the source directory
  for (fs::directory_iterator file(source); file != fs::directory_iterator();
       ++file) {
    try {
      fs::path current(file->path());
      if (fs::is_directory(current)) {
        // Found directory: Recursion
        if (!copyDir(current, destination / current.filename())) {
          return false;
        }
      } else {
        // Found file: Copy
        fs::copy_file(current, destination / current.filename());
      }
    } catch (fs::filesystem_error const &e) {
      std::cerr << e.what() << std::endl;
    }
  }
  return true;
}

// removes existing backup folder and creates a new backup
void backUpDir(const std::string &_cur_dir, const std::string &_bak_dir) {
  boost::filesystem::path cur_dir(_cur_dir);
  boost::filesystem::path bak_dir(_bak_dir);

  std::cout << "Backing up " << cur_dir << std::endl;
  if (boost::filesystem::is_directory(bak_dir)) {
    std::cout << "Removing existing " << bak_dir << std::endl;
    int nRem = boost::filesystem::remove_all(bak_dir);
    std::cout << "Removed " << nRem << " item(s)." << std::endl;
  }

  // recursive copying
  if (!copyDir(cur_dir, bak_dir)) {
    std::cerr << "Problem while copying directory " << _cur_dir << std::endl;
    throw;
  }
}

// remove files inside a directory
void rmDirContent(const boost::filesystem::path &_dir, bool _skipDir = true) {
  for (boost::filesystem::directory_iterator end_dir_it, it(_dir);
       it != end_dir_it; ++it) {
    if (boost::filesystem::is_directory(it->path()) && _skipDir) continue;
    boost::filesystem::remove_all(it->path());
  }
}

// remove recursive all content that contain given name
void rmFileName(const boost::filesystem::path &_dir, const std::string &_fn) {
  // iterate and delete matches
  for (boost::filesystem::directory_iterator end_dir_it, it(_dir);
       it != end_dir_it; ++it) {
    std::string fext(it->path().filename().string());
    if (!expContains(fext, _fn)) continue;
    std::cout << "Deleting " << it->path() << std::endl;
    boost::filesystem::remove_all(it->path());
  }
}

// check if name contains _exp
bool expContains(const std::string &_name, const std::string &_exp) {
  return _name.find(_exp) != std::string::npos;
}

// recursive copy all content that contain given name
void cpFileName(const boost::filesystem::path &_src_dir,
                const boost::filesystem::path &_des_dir, const std::string &_fn,
                const std::string &_starts) {
  // iterate and delete matches
  for (boost::filesystem::directory_iterator end_dir_it, it(_src_dir);
       it != end_dir_it; ++it) {
    bool overide = false;
    if (expContains(_starts, "*")) overide = true;

    std::string fext(it->path().filename().string());
    if (!expContains(fext, _fn)) continue;
    {
      if (!overide)
        if (fext.compare(0, _starts.size(), _starts)) continue;

      boost::filesystem::path des_path(_des_dir.string() + "/" +
                                       it->path().filename().string());
      std::cout << "Copying " << it->path() << " to " << des_path << std::endl;
      boost::filesystem::copy_file(
          it->path(), des_path,
          boost::filesystem::copy_option::overwrite_if_exists);
    }
  }
}

}  // namespace DRV
}  // namespace NEM

#endif  // HAVE_SIMMETRIX
