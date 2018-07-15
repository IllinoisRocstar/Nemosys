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

RemeshDriver::RemeshDriver(const std::vector<std::string>& _fluidNames,
                           const std::vector<std::string>& _ifluidniNames,
                           const std::vector<std::string>& _ifluidnbNames,
                           const std::vector<std::string>& _ifluidbNames,
                           const std::vector<std::string>& _burnNames,
                           const std::vector<std::string>& _iBurnNames,
                           const json& remeshjson, int _numPartitions)
  : fluidNames(_fluidNames), ifluidniNames(_ifluidniNames), ifluidnbNames(_ifluidnbNames),
    ifluidbNames(_ifluidbNames), burnNames(_burnNames), iBurnNames(_iBurnNames),
    numPartitions(_numPartitions)
{
  // stitch fluid files
  stitchCGNS(fluidNames,0);
	// stitch burn files
	stitchCGNS(burnNames,0);
	// stitch iburn files
	stitchCGNS(iBurnNames,1);
  // stitch ifluid_ni files
  stitchCGNS(ifluidniNames,1);
  // stitch ifluid_b files
  stitchCGNS(ifluidbNames,1);
  // stitch ifluid_ng files
  stitchCGNS(ifluidnbNames,1);
  // get stitched meshes from stitcher vector
  for (int i = 0; i < stitchers.size(); ++i)
  {
    //cgObjs.push_back(stitchers[i]->getStitchedCGNS());
    mbObjs.push_back(stitchers[i]->getStitchedMB());
    mbObjs[i]->setContBool(0);
  }
  // creates remeshedVol and remeshedSurf
  remesh(remeshjson);
  // creates stitchedSurf
  if (mbObjs.size() > 1)
  {
    // stitches b, ni and nb surfaces into one stitchedSurf
    stitchSurfaces();
    // do not smooth for cell data transfer
    stitchedSurf->setContBool(0);
    // transfer patch number to remeshed surface
    std::vector<std::string> patchno(1,"patchNo");
    stitchedSurf->transfer(remeshedSurf.get(), "Consistent Interpolation", patchno, 1);
    stitchedSurf->write();
  }
  remeshedSurf->write();
  std::unique_ptr<RocPartCommGenDriver> rocprepdrvr 
    = std::unique_ptr<RocPartCommGenDriver>
        (new RocPartCommGenDriver(this->remeshedVol, this->remeshedSurf, numPartitions));
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
  std::cout << "Extracting surface mesh #############################################\n";
  std::unique_ptr<meshBase> surf 
    = std::unique_ptr<meshBase>
        (meshBase::Create(mbObjs[0]->extractSurface(), "extractedSurface.vtp")); 
  // remeshing with engine specified in input
  surf->write("extractedSurface.stl"); 
  mshgendrvr 
    = std::unique_ptr<MeshGenDriver>
        (MeshGenDriver::readJSON("extractedSurface.stl", "remeshedVol.vtu", remeshjson));
  // get the remeshed volume from the mesh generator
  remeshedVol = mshgendrvr->getNewMesh();
  // extract the remeshed volumes surface  
  remeshedSurf = meshBase::CreateShared(remeshedVol->extractSurface(), "remeshedSurf.vtp");
}

void RemeshDriver::stitchSurfaces()
{
  std::cout << mbObjs.size() << std::endl;
  // stitch b, ni and nb surfaces
  std::vector<std::shared_ptr<meshBase>> surfs;
  surfs.insert(surfs.begin(), mbObjs.begin()+1,mbObjs.end());
  stitchedSurf = meshBase::stitchMB(surfs);
  stitchedSurf->setFileName("stitchedSurf.vtu");
}

RemeshDriver* RemeshDriver::readJSON(json inputjson)
{
  std::string case_dir = inputjson["RocFluMP Case Directory"].as<std::string>();
  std::string burn_dir = inputjson["RocFluMP Case Directory"].as<std::string>();
  std::string base_t = inputjson["Base Time Step"].as<std::string>();
  std::vector<std::string> fluNames(getCgFNames(case_dir, "fluid", base_t));
  std::vector<std::string> ifluniNames(getCgFNames(case_dir, "ifluid_ni", base_t));
  std::vector<std::string> iflunbNames(getCgFNames(case_dir, "ifluid_nb", base_t));
  std::vector<std::string> iflubNames(getCgFNames(case_dir, "ifluid_b", base_t));
  std::vector<std::string> burnNames(getCgFNames(burn_dir, "burn", base_t));
	std::vector<std::string> iBurnNames(getCgFNames(burn_dir, "iburn", base_t));
  json remeshjson = inputjson["Remeshing Options"];
  int numPartitions = inputjson["Number Of New Partitions"].as<int>();
  RemeshDriver* remeshdrvobj = new RemeshDriver(fluNames, ifluniNames, iflunbNames, 
                                                iflubNames, burnNames, iBurnNames,
                                                remeshjson, numPartitions);
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
