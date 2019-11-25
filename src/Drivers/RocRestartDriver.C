#include <RocRestartDriver.H>
#include <meshStitcher.H>
#include <rocstarCgns.H>
#include <AuxiliaryFunctions.H>

#include "TransferDriver.H"

RocRestartDriver::RocRestartDriver(const std::vector<std::string>& _fluidNamesRm,
                                   const std::vector<std::string>& _ifluidniNamesRm,
                                   const std::vector<std::string>& _ifluidnbNamesRm,
                                   const std::vector<std::string>& _ifluidbNamesRm,
                                   const std::vector<std::string>& _fluidNamesLts,
                                   const std::vector<std::string>& _ifluidniNamesLts,  
                                   const std::vector<std::string>& _ifluidnbNamesLts,   
                                   const std::vector<std::string>& _ifluidbNamesLts,
                                   const std::vector<std::string>& _burnNamesRm,
                                   const std::vector<std::string>& _iBurnNamesRm, 
                                   const std::vector<std::string>& _burnNamesLts, 
                                   const std::vector<std::string>& _iBurnNamesLts)
  : fluidNamesRm(_fluidNamesRm), ifluidniNamesRm(_ifluidniNamesRm),
    ifluidnbNamesRm(_ifluidnbNamesRm), ifluidbNamesRm(_ifluidbNamesRm),
    fluidNamesLts(_fluidNamesLts), ifluidniNamesLts(_ifluidniNamesLts),
    ifluidnbNamesLts(_ifluidnbNamesLts), ifluidbNamesLts(_ifluidbNamesLts),
	burnNamesRm(_burnNamesRm), iBurnNamesRm(_iBurnNamesRm), 
	burnNamesLts(_iBurnNamesLts), iBurnNamesLts(_iBurnNamesLts)
{
  //---- stitch files from last time step
  
  // fluid
  stitchCGNS(fluidNamesLts,0);
  // burn
  stitchCGNS(burnNamesLts,0);
  // iburn
  stitchCGNS(iBurnNamesLts,1);
  // ifluid_ni 
  stitchCGNS(ifluidniNamesLts,1);
  // ifluid_b 
  stitchCGNS(ifluidbNamesLts,1);
  // ifluid_ng 
  stitchCGNS(ifluidnbNamesLts,1);
  // get stitched meshes from stitcher vector
  for (int i = 0; i < stitchers.size(); ++i)
  {
    mbObjs.push_back(stitchers[i]->getStitchedMB());
    mbObjs[i]->setContBool(0);
  }
  // stitch b, ni and nb surfaces to create stitchedSurf
  stitchSurfaces(); 

  //---- load remeshed cgns partitions into vecs and populate converted MB vecs
  loadPartCgMb();

  //---- transfer solution fields from stitchedMBs to each MB partition 
  //     and push into corresponding open cgns file
  transferStitchedToPartCg("Consistent Interpolation");

  std::cout << "RocRestartDriver created" << std::endl; 
} 


RocRestartDriver::~RocRestartDriver()
{
  std::cout << "RocRestartDriver destroyed" << std::endl;
}

void RocRestartDriver::stitchCGNS(const std::vector<std::string>& fnames, bool surf)
{
  if (fnames.size())
  {
    stitchers.push_back(std::unique_ptr<meshStitcher>(new meshStitcher(fnames, surf)));
  }
}

void RocRestartDriver::stitchSurfaces()
{
  // stitch b, ni and nb surfaces
  std::vector<std::shared_ptr<meshBase>> surfs;
  surfs.insert(surfs.begin(), mbObjs.begin()+3,mbObjs.end());
  stitchedSurf = meshBase::stitchMB(surfs);
  stitchedSurf->setContBool(0);
  stitchedSurf->setFileName("stitchedSurf.vtu");
}

void RocRestartDriver::loadPartCgMb()
{
  // fluid
  if (fluidNamesRm.size())
  {
    cgVtPair fluidPair(loadCGNS(fluidNamesRm,0));
    fluidRmCg = fluidPair.first;
    fluidRmMb = fluidPair.second;
  }
	// burn
	if (burnNamesRm.size())
	{
		cgVtPair burnPair(loadCGNS(burnNamesRm,0));
		burnRmCg = burnPair.first;
		burnRmMb = burnPair.second;
	}	
	// iburn
	if (iBurnNamesRm.size())
	{
		cgVtPair iBurnPair(loadCGNS(iBurnNamesRm,1));
		iBurnRmCg = iBurnPair.first;
		iBurnRmMb = iBurnPair.second;
	}
  // ifluid_ni
  if (ifluidniNamesRm.size())
  {
    cgVtPair ifluidNiPair(loadCGNS(ifluidniNamesRm,1));
    ifluidNiRmCg = ifluidNiPair.first;
    ifluidNiRmMb = ifluidNiPair.second;
  }
  // ifluid_nb
  if (ifluidnbNamesRm.size())
  {
    cgVtPair ifluidNbPair(loadCGNS(ifluidnbNamesRm,1));
    ifluidNbRmCg = ifluidNbPair.first;
    ifluidNbRmMb = ifluidNbPair.second;
  }
  // ifluid_b
  if (ifluidbNamesRm.size())
  {
    cgVtPair ifluidBPair(loadCGNS(ifluidbNamesRm,1));
    ifluidBRmCg = ifluidBPair.first;
    ifluidBRmMb = ifluidBPair.second;
  }
}

void RocRestartDriver::transferStitchedToPartCg(const std::string& transferType)
{
  // fluid
  for (int i = 0; i < fluidRmMb.size(); ++i)
  {
    // mbObjs[0]->transfer(fluidRmMb[i].get(), transferType); 
    auto transfer = TransferDriver::CreateTransferObject(mbObjs[0], fluidRmMb[i].get(), transferType);
    transfer->run(mbObjs[0]->getNewArrayNames());

    fluidRmCg[i]->overwriteSolData(fluidRmMb[i].get());
  }
	// burn
	for (int i = 0; i < burnRmMb.size(); ++i)
	{
		// mbObjs[1]->transfer(burnRmMb[i].get(), transferType);
    auto transfer = TransferDriver::CreateTransferObject(mbObjs[1], burnRmMb[i].get(), transferType);
    transfer->run(mbObjs[1]->getNewArrayNames());

		burnRmCg[i]->overwriteSolData(burnRmMb[i].get());
	}
	// iburn
	for (int i = 0; i < iBurnRmMb.size(); ++i)
	{
		// mbObjs[2]->transfer(iBurnRmMb[i].get(), transferType);
    auto transfer = TransferDriver::CreateTransferObject(mbObjs[2], iBurnRmMb[i].get(), transferType);
    transfer->run(mbObjs[2]->getNewArrayNames());

		iBurnRmCg[i]->overwriteSolData(iBurnRmMb[i].get());
	}
  // ifluid_ni
  for (int i = 0; i < ifluidNiRmMb.size(); ++i)
  {
    // stitchedSurf->transfer(ifluidNiRmMb[i].get(), transferType);
    auto transfer = TransferDriver::CreateTransferObject(stitchedSurf, ifluidNiRmMb[i].get(), transferType);
    transfer->run(stichedSurf->getNewArrayNames());

    ifluidNiRmCg[i]->overwriteSolData(ifluidNiRmMb[i].get());
  }
  // ifluid_nb
  for (int i = 0; i < ifluidNbRmMb.size(); ++i)
  {
    // stitchedSurf->transfer(ifluidNbRmMb[i].get(), transferType);
    auto transfer = TransferDriver::CreateTransferObject(stitchedSurf, ifluidNbRmMb[i].get(), transferType);
    transfer->run(stichedSurf->getNewArrayNames());

    ifluidNbRmCg[i]->overwriteSolData(ifluidNbRmMb[i].get());
  }
  // ifluid_b
  for (int i = 0; i < ifluidBRmMb.size(); ++i)
  {
    // stitchedSurf->transfer(ifluidBRmMb[i].get(), transferType);
    auto transfer = TransferDriver::CreateTransferObject(stitchedSurf, ifluidBRmMb[i].get(), transferType);
    transfer->run(stichedSurf->getNewArrayNames());

    ifluidBRmCg[i]->overwriteSolData(ifluidBRmMb[i].get());
  }
}


cgVtPair RocRestartDriver::loadCGNS(const std::vector<std::string>& fnames, bool surf)
{
  std::vector<std::shared_ptr<cgnsAnalyzer>> cgObjs;
  std::vector<std::shared_ptr<meshBase>> mbobjs(fnames.size());
  for (int i = 0; i < fnames.size(); ++i)
  {
    std::shared_ptr<cgnsAnalyzer> cgObj;
    if (surf)
    {
      cgnsAnalyzer* _cgObj = new rocstarCgns(fnames[i]);
      cgObj.reset(_cgObj);
    }
    else
    {
      cgnsAnalyzer* _cgObj = new cgnsAnalyzer(fnames[i]);
      cgObj.reset(_cgObj);
    }
    cgObj->loadGrid(1);
    cgObjs.push_back(cgObj);
    std::size_t pos = fnames[i].find_last_of("/");
    std::string vtkname = fnames[i].substr(pos+1);
    vtkname = nemAux::trim_fname(vtkname, ".vtu");
    mbobjs[i] = meshBase::CreateShared(cgObj->getVTKMesh(),vtkname);      
  }
  return std::make_pair(cgObjs, mbobjs);
}

RocRestartDriver* RocRestartDriver::readJSON(json inputjson)
{
  std::string prepDir = inputjson["Rocprep Remesh Directory"].as<std::string>();
	std::string prepBurnDir = inputjson["Rocprep Remesh Burn Directory"].as<std::string>();
  std::string lastDir = inputjson["Rocout Last TS Directory"].as<std::string>();
	std::string lastBurnDir = inputjson["Rocout Last TS Burn Directory"].as<std::string>();
  std::string base_t = inputjson["Rocout Last TS Base"].as<std::string>();
  std::string base_tRm = inputjson["Rocout Remesh TS Base"].as<std::string>();

  std::vector<std::string> fluNamesRm(getCgFNames(prepDir, "fluid", base_tRm));
  std::vector<std::string> burnNamesRm(getCgFNames(prepBurnDir, "burn", base_tRm));
	std::vector<std::string> iBurnNamesRm(getCgFNames(prepBurnDir, "iburn", base_tRm));
	std::vector<std::string> ifluniNamesRm(getCgFNames(prepDir, "ifluid_ni", base_tRm));
  std::vector<std::string> iflunbNamesRm(getCgFNames(prepDir, "ifluid_nb", base_tRm));
  std::vector<std::string> iflubNamesRm(getCgFNames(prepDir, "ifluid_b", base_tRm));
  std::vector<std::string> fluNamesLts(getCgFNames(lastDir, "fluid", base_t));
  std::vector<std::string> burnNamesLts(getCgFNames(lastBurnDir, "burn", base_t));
	std::vector<std::string> iBurnNamesLts(getCgFNames(lastBurnDir, "iburn", base_t));
  std::vector<std::string> ifluniNamesLts(getCgFNames(lastDir, "ifluid_ni", base_t));
  std::vector<std::string> iflunbNamesLts(getCgFNames(lastDir, "ifluid_nb", base_t));
  std::vector<std::string> iflubNamesLts(getCgFNames(lastDir, "ifluid_b", base_t));
  
  RocRestartDriver* restartDriver 
    = new RocRestartDriver(fluNamesRm, ifluniNamesRm, iflunbNamesRm, iflubNamesRm,
                           fluNamesLts, ifluniNamesLts, iflunbNamesLts, iflubNamesLts,
													 burnNamesRm, iBurnNamesRm, burnNamesLts, iBurnNamesLts);
  return restartDriver;
}

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
}
