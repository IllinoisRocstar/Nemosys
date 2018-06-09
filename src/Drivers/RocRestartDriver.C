#include <RocRestartDriver.H>
#include <meshStitcher.H>
#include <rocstarCgns.H>

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

void RocRestartDriver::deleteCgMb(std::vector<cgnsAnalyzer*>& cgObjs, std::vector<meshBase*>& cgMbObjs)
{
  for (int i = 0; i < cgObjs.size(); ++i)
  {
    if (cgObjs[i])
    {
      delete cgObjs[i];
      delete cgMbObjs[i];
      cgObjs[i] = nullptr;
      cgMbObjs[i] = nullptr;
    }
  }
}

RocRestartDriver::~RocRestartDriver()
{
  deleteCgMb(fluidRmCg, fluidRmMb);
  deleteCgMb(ifluidNiRmCg, ifluidNiRmMb);
  deleteCgMb(ifluidBRmCg, ifluidBRmMb);
  deleteCgMb(ifluidNbRmCg, ifluidNbRmMb);
	deleteCgMb(burnRmCg, burnRmMb);
	deleteCgMb(iBurnRmCg, iBurnRmMb); 
 
  for (int i = 0; i < stitchers.size(); ++i)
  {
    if (stitchers[i])
    { 
      delete stitchers[i]; // deletes mb and cg objs as well
      stitchers[i] = nullptr;
      mbObjs[i] = nullptr;
    }
  }

  if (stitchedSurf)
  {
    delete stitchedSurf;
    stitchedSurf = nullptr;
  } 
  std::cout << "RocRestartDriver destroyed" << std::endl;
}

void RocRestartDriver::stitchCGNS(const std::vector<std::string>& fnames, bool surf)
{
  if (fnames.size())
  {
    stitchers.push_back(new meshStitcher(fnames, surf));
  }
}

void RocRestartDriver::stitchSurfaces()
{
  // stitch b, ni and nb surfaces
  std::vector<meshBase*> surfs;
  surfs.insert(surfs.begin(), mbObjs.begin()+3,mbObjs.end());
  stitchedSurf = meshBase::stitchMB(surfs);
  stitchedSurf->setContBool(0);
  stitchedSurf->setFileName("stitchedSurf.vtu");
}

void RocRestartDriver::loadPartCgMb()
{
  // return type from loadCGNS 
  typedef std::pair<std::vector<cgnsAnalyzer*>, std::vector<meshBase*>> cgVtPair; 
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
    mbObjs[0]->transfer(fluidRmMb[i], transferType); 
    fluidRmCg[i]->overwriteSolData(fluidRmMb[i]);
  }
	// burn
	for (int i = 0; i < burnRmMb.size(); ++i)
	{
		mbObjs[1]->transfer(burnRmMb[i], transferType);
		burnRmCg[i]->overwriteSolData(burnRmMb[i]);
	}
	// iburn
	for (int i = 0; i < iBurnRmMb.size(); ++i)
	{
		mbObjs[2]->transfer(iBurnRmMb[i], transferType);
		iBurnRmCg[i]->overwriteSolData(iBurnRmMb[i]);
	}
  // ifluid_ni
  for (int i = 0; i < ifluidNiRmMb.size(); ++i)
  {
    stitchedSurf->transfer(ifluidNiRmMb[i], transferType);
    ifluidNiRmCg[i]->overwriteSolData(ifluidNiRmMb[i]);
  }
  // ifluid_nb
  for (int i = 0; i < ifluidNbRmMb.size(); ++i)
  {
    stitchedSurf->transfer(ifluidNbRmMb[i], transferType);
    ifluidNbRmCg[i]->overwriteSolData(ifluidNbRmMb[i]);
  }
  // ifluid_b
  for (int i = 0; i < ifluidBRmMb.size(); ++i)
  {
    stitchedSurf->transfer(ifluidBRmMb[i], transferType);
    ifluidBRmCg[i]->overwriteSolData(ifluidBRmMb[i]);
  }
}

std::pair< std::vector<cgnsAnalyzer*>, std::vector<meshBase*> >
RocRestartDriver::loadCGNS(const std::vector<std::string>& fnames, bool surf)
{
  std::vector<cgnsAnalyzer*> cgObjs;
  std::vector<meshBase*> mbobjs(fnames.size());
  for (int i = 0; i < fnames.size(); ++i)
  {
    cgnsAnalyzer* cgObj;
    if (surf)
    {
      cgObj = new rocstarCgns(fnames[i]);
    }
    else
    {
      cgObj = new cgnsAnalyzer(fnames[i]);
    }
    cgObj->loadGrid(1);
    cgObjs.push_back(cgObj);
    std::size_t pos = fnames[i].find_last_of("/");
    std::string vtkname = fnames[i].substr(pos+1);
    vtkname = trim_fname(vtkname, ".vtu");
    mbobjs[i] = meshBase::Create(cgObj->getVTKMesh(),vtkname);      
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
