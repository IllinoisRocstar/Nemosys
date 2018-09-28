// nemosys
#include <RemeshDriver.H>
#include <meshStitcher.H>
#include <rocstarCgns.H>
#include <cgnsWriter.H> // not used
#include <meshPartitioner.H> // not used
#include <MeshGenDriver.H>
#include <RocPartCommGenDriver.H>
#include <AuxiliaryFunctions.H>

//vtk
#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGenericCell.h>
#include <vtkCell.h>
#include <vtkCellData.h>

// file system
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

void backUpDir(std::string _cur_dir, std::string _bak_dir);
void rmDirContent(boost::filesystem::path const & _dir, bool _skipDir);
bool expContains(std::string _name, std::string _exp);
void rmFileName(boost::filesystem::path const & _dir, std::string _fn);
void cpFileName(boost::filesystem::path const & _src_dir, boost::filesystem::path const & _des_dir, std::string _fn, std::string _starts);


RemeshDriver::RemeshDriver(const std::vector<std::string>& _fluidNames,
                           const std::vector<std::string>& _ifluidniNames,
                           const std::vector<std::string>& _ifluidnbNames,
                           const std::vector<std::string>& _ifluidbNames,
                           const std::vector<std::string>& _burnNames,
                           const std::vector<std::string>& _iBurnNames,
                           const json& remeshjson, int _numPartitions,
                           const std::string& base_t, int writeIntermediateFiles,
                           double searchTolerance, const std::string& caseName,
                           const std::string _prefix_path)
  : fluidNames(_fluidNames), ifluidniNames(_ifluidniNames), ifluidnbNames(_ifluidnbNames),
    ifluidbNames(_ifluidbNames), burnNames(_burnNames), iBurnNames(_iBurnNames),
    numPartitions(_numPartitions)
{
  // setting path
  prefixPath = (_prefix_path.length()!=0) ? _prefix_path+"/":"./";

  // preparing input for the actural operator class
  // stitch fluid files
  this->stitchCGNS(fluidNames,0);
  // stitch burn files
  this->stitchCGNS(burnNames,0);
  // stitch iburn files
  this->stitchCGNS(iBurnNames,1);
  // stitch ifluid_ni files
  this->stitchCGNS(ifluidniNames,1);
  // stitch ifluid_b files
  this->stitchCGNS(ifluidbNames,1);
  // stitch ifluid_nb files
  this->stitchCGNS(ifluidnbNames,1);
  // get stitched meshes from stitcher vector
  for (int i = 0; i < this->stitchers.size(); ++i)
  {
    //cgObjs.push_back(stitchers[i]->getStitchedCGNS());
    this->mbObjs.push_back(stitchers[i]->getStitchedMB());
    this->mbObjs[i]->setContBool(0);
  }
  if (writeIntermediateFiles) this->mbObjs[0]->write(prefixPath+"stitchedVol.vtu");
  // creates remeshedVol and remeshedSurf
  this->remesh(remeshjson);
  // creates stitchedSurf
  if (this->mbObjs.size() > 1)
  {
    // cgns files are already stitched together into meshbase obj
    // need to transfer this data onto the surface mesh

    // stitches b, ni and nb surfaces into one stitchedSurf
    // also stitches iburn files
    this->stitchSurfaces();
    // do not smooth for cell data transfer
    this->stitchedSurf->setContBool(0);
    this->stitchedBurnSurf->setContBool(0);
    // transfer patch number to remeshed surface
    std::vector<std::string> transferPaneData{"patchNo", "bcflag", "cnstr_type"};
    this->stitchedSurf->transfer(remeshedSurf.get(), "Consistent Interpolation", transferPaneData, 1);
    if (writeIntermediateFiles) this->stitchedSurf->write();
    if (writeIntermediateFiles) this->stitchedBurnSurf->write();
  }
  if (writeIntermediateFiles) this->remeshedSurf->write();
  bool withC2CTransSmooth = 
      (remeshjson.has_key("C2CTransSmooth") ? remeshjson["C2CTransSmooth"].as<bool>() : false);

  // instantiating and executing actual operator class
  std::unique_ptr<RocPartCommGenDriver> rocprepdrvr 
    = std::unique_ptr<RocPartCommGenDriver>
        (
         new RocPartCommGenDriver(this->remeshedVol, this->remeshedSurf, 
                                  this->mbObjs[0], this->stitchedSurf,
                                  this->stitchedBurnSurf,
                                  numPartitions, base_t, writeIntermediateFiles, 
                                  searchTolerance, caseName,
                                  withC2CTransSmooth, prefixPath) 
        );
  std::cout << "RemeshDriver created" << std::endl;
}
 
RemeshDriver::~RemeshDriver()
{
  std::cout << "RemeshDriver destroyed" << std::endl;
}


void RemeshDriver::stitchCGNS(const std::vector<std::string>& fnames, bool surf)
{
  if (fnames.size())
  {
    stitchers.push_back(std::unique_ptr<meshStitcher>(new meshStitcher(fnames, surf)));
  }
}

void RemeshDriver::remesh(const json& remeshjson)
{
  std::cout << "Extracting surface mesh\n";
  std::unique_ptr<meshBase> surf = std::unique_ptr<meshBase>
        (meshBase::Create(mbObjs[0]->extractSurface(), prefixPath+"extractedSurface.vtp")); 
  // remeshing with engine specified in input
  surf->write(prefixPath+"extractedSurface.stl"); 
  mshgendrvr = std::unique_ptr<MeshGenDriver>
        (MeshGenDriver::readJSON(prefixPath+"extractedSurface.stl", prefixPath+"remeshedVol.vtu", remeshjson));
  // get the remeshed volume from the mesh generator
  remeshedVol = mshgendrvr->getNewMesh();
  // extract the remeshed volumes surface  
  remeshedSurf = meshBase::CreateShared(remeshedVol->extractSurface(), prefixPath+"remeshedSurf.vtp");
}

void RemeshDriver::stitchSurfaces()
{
  std::cout << mbObjs.size() << std::endl;
  // stitch b, ni and nb surfaces
  std::vector<std::shared_ptr<meshBase>> surfs;
  surfs.insert(surfs.begin(), mbObjs.begin()+3,mbObjs.end());
  stitchedSurf = meshBase::stitchMB(surfs);
  stitchedSurf->setFileName(prefixPath+"stitchedSurf.vtu");

  std::vector<std::shared_ptr<meshBase>> burnSurfs;
  burnSurfs.insert(burnSurfs.begin(), mbObjs.begin()+2, mbObjs.begin()+3);
  std::cout << "stitching burning surfaces" << std::endl;
  stitchedBurnSurf = meshBase::stitchMB(burnSurfs);
  stitchedBurnSurf->setFileName(prefixPath+"stitchedBurnSurf.vtu");
}

RemeshDriver* RemeshDriver::readJSON(json inputjson)
{
  std::string case_name = inputjson["Rocstar Case Name"].as<std::string>();
  std::string case_dir = inputjson["RocFluMP Case Directory"].as<std::string>();
  std::string burn_dir = inputjson["RocBurnAPN Case Directory"].as<std::string>();
  std::string base_t = inputjson["Base Time Step"].as<std::string>();
  std::vector<std::string> fluNames(getCgFNames(case_dir, "fluid", base_t));
  std::vector<std::string> ifluniNames(getCgFNames(case_dir, "ifluid_ni", base_t));
  std::vector<std::string> iflunbNames(getCgFNames(case_dir, "ifluid_nb", base_t));
  std::vector<std::string> iflubNames(getCgFNames(case_dir, "ifluid_b", base_t));
  std::vector<std::string> burnNames(getCgFNames(burn_dir, "burn", base_t));
  std::vector<std::string> iBurnNames(getCgFNames(burn_dir, "iburn", base_t));
  if (fluNames.size()==0 || burnNames.size()==0)
  {
    std::cerr << "Error finding Rocstar output files. Make sure folder locations are properly set in the input file and all folders are non-empty.\n";
    throw;
  }
  json remeshjson = inputjson["Remeshing Options"];
  int numPartitions = inputjson["Number Of New Partitions"].as<int>();
  int writeIntermediateFiles = 0;
  double searchTolerance = 1e-9;
  if (inputjson.has_key("Search Tolerance"))
  {
    searchTolerance = inputjson["Search Tolerance"].as<double>();
  }
  if (inputjson.has_key("Write Intermediate Files"))
  {
    writeIntermediateFiles = inputjson["Write Intermediate Files"].as<int>();
  }
  
  // set temporary dir
  const char dir_path[] = "/remesh_tmp";

  // preforming prep actions
  boost::filesystem::path cwd_dir;
  cwd_dir = boost::filesystem::current_path();
  boost::filesystem::path rmsh_dir(cwd_dir);
  rmsh_dir += boost::filesystem::path(dir_path);  
  std::cout << "Current working directory " << cwd_dir << "\n";

  // make a tmp directory or clean if already there
  if (boost::filesystem::is_directory(rmsh_dir))
  {
      std::cout << "Removing existing " << rmsh_dir << std::endl;
      int nRem = boost::filesystem::remove_all(rmsh_dir);
      std::cout << "Removed " << nRem << " item(s).\n";
      if (boost::filesystem::create_directory(rmsh_dir))
          std::cout << "Created " << rmsh_dir << std::endl;
      else
      {
          std::cerr << "Problem creating directory " << rmsh_dir << "\n";
          throw;
      }
  }
  else if (boost::filesystem::create_directory(rmsh_dir)) 
  {
      std::cout << "Created " << rmsh_dir << std::endl;
  }

  // backing up all old cgns files for rocflu (TODO: expensive remove)
  backUpDir(case_dir, rmsh_dir.string()+"/_flu_bak");
  backUpDir(burn_dir, rmsh_dir.string()+"/_burn_bak");
  

  // remeshing, repartitioning and solution transfer actions
  RemeshDriver* remeshdrvobj = new RemeshDriver(fluNames, ifluniNames, iflunbNames, 
                                                iflubNames, burnNames, iBurnNames,
                                                remeshjson, numPartitions, base_t,
                                                writeIntermediateFiles, searchTolerance,
                                                case_name, rmsh_dir.string());
  // post action
  // remove old cgns files (TODO: copy somewhere)
  // remove old cgns and files
  rmFileName(case_dir, ".cgns");
  rmFileName(case_dir, ".txt");
  rmFileName(burn_dir, ".cgns");
  rmFileName(burn_dir, ".txt");

  // remove flu-modin files
  boost::filesystem::path flumodin_dir(case_dir+"/../Modin/");  
  rmFileName(flumodin_dir, ".map");
  rmFileName(flumodin_dir, ".dim");
  rmFileName(flumodin_dir, ".com");
  rmFileName(flumodin_dir, ".cmp");
  rmFileName(flumodin_dir, ".rnm");

  // copy over old cgns and txt files
  cpFileName(rmsh_dir.string()+"/_flu_bak/", case_dir, "fluid", "fluid");
  cpFileName(rmsh_dir.string()+"/_flu_bak/", case_dir, "ifluid", "ifluid");
  cpFileName(rmsh_dir.string()+"/_burn_bak/", burn_dir, "burn", "burn");
  cpFileName(rmsh_dir.string()+"/_burn_bak/", burn_dir, "iburn", "iburn");


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

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
}

// file transfer helpers
// recursive copy, will create distenation folder (make sure it is removed first)
bool copyDir(boost::filesystem::path const & source, boost::filesystem::path const & destination
)
{
    namespace fs = boost::filesystem;
    try
    {
        // Check whether the function call is valid
        if(
            !fs::exists(source) ||
            !fs::is_directory(source)
        )
        {
            std::cerr << "Source directory " << source.string()
                << " does not exist or is not a directory." << '\n'
            ;
            return false;
        }
        if(fs::exists(destination))
        {
            std::cerr << "Destination directory " << destination.string()
                << " already exists." << '\n'
            ;
            return false;
        }
        // Create the destination directory
        if(!fs::create_directory(destination))
        {
            std::cerr << "Unable to create destination directory"
                << destination.string() << '\n'
            ;
            return false;
        }
    }
    catch(fs::filesystem_error const & e)
    {
        std::cerr << e.what() << '\n';
        return false;
    }
    // Iterate through the source directory
    for(
        fs::directory_iterator file(source);
        file != fs::directory_iterator(); ++file
    )
    {
        try
        {
            fs::path current(file->path());
            if(fs::is_directory(current))
            {
                // Found directory: Recursion
                if(
                    !copyDir(
                        current,
                        destination / current.filename()
                    )
                )
                {
                    return false;
                }
            }
            else
            {
                // Found file: Copy
                fs::copy_file(
                    current,
                    destination / current.filename()
                );
            }
        }
        catch(fs::filesystem_error const & e)
        {
            std:: cerr << e.what() << '\n';
        }
    }
    return true;
}

// removes existing backup folder and creates a new backup
void backUpDir(std::string _cur_dir, std::string _bak_dir)
{
  boost::filesystem::path cur_dir(_cur_dir);
  boost::filesystem::path bak_dir(_bak_dir);

  std::cout << "Backing up " << cur_dir << "\n";
  if (boost::filesystem::is_directory(bak_dir))
  {
      std::cout << "Removing existing " << bak_dir << std::endl;
      int nRem = boost::filesystem::remove_all(bak_dir);
      std::cout << "Removed " << nRem << " item(s).\n";
  }

  // recursive copying
  if (!copyDir(cur_dir, bak_dir))
  {
      std::cerr << "Problem while copying directory " << _cur_dir << "\n";
      throw;
  }
}

// remove files inside a directory
void rmDirContent(boost::filesystem::path const & _dir, bool _skipDir = true)
{
  for (boost::filesystem::directory_iterator end_dir_it, it(_dir); it!=end_dir_it; ++it) 
  {
      if (boost::filesystem::is_directory(it->path()) && _skipDir)
          continue;
      boost::filesystem::remove_all(it->path());
  }
}

// remove recursive all content that contain given name
void rmFileName(boost::filesystem::path const & _dir, std::string _fn)
{
  // iterate and delete matches
  for (boost::filesystem::directory_iterator end_dir_it, it(_dir); it!=end_dir_it; ++it) 
  {
      std::string fext(it->path().filename().string());
      if ( !expContains(fext, _fn) ) continue;
      std::cout << "Deleting " << it->path() << std::endl;
      boost::filesystem::remove_all(it->path());
  }
}

// check if name conatins _exp
bool expContains(std::string _name, std::string _exp)
{
    return (_name.find(_exp) != std::string::npos);
}

// recursive copy all content that contain given name
void cpFileName(boost::filesystem::path const & _src_dir, boost::filesystem::path const & _des_dir, std::string _fn, std::string _starts)
{
  // iterate and delete matches
  for (boost::filesystem::directory_iterator end_dir_it, it(_src_dir); it!=end_dir_it; ++it) 
  {
      bool overide = false;
      if (expContains(_starts,"*"))
          overide = true;

      std::string fext(it->path().filename().string());
      if ( !expContains(fext, _fn) ) continue;
      {
        if (!overide)  
            if ( fext.compare(0, _starts.size(), _starts) )
                continue;

        boost::filesystem::path des_path(_des_dir.string() + "/"+it->path().filename().string() );
        std::cout << "Copying " << it->path() 
                  << " to " << des_path << std::endl;
        boost::filesystem::copy_file(it->path(), des_path, boost::filesystem::copy_option::overwrite_if_exists);
      }
  }
}

