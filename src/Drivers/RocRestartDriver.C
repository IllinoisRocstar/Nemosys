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
                                   const std::vector<std::string>& _ifluidbNamesLts)
  : fluidNamesRm(_fluidNamesRm), ifluidniNamesRm(_ifluidniNamesRm),
    ifluidnbNamesRm(_ifluidnbNamesRm), ifluidbNamesRm(_ifluidbNamesRm),
    fluidNamesLts(_fluidNamesLts), ifluidniNamesLts(_ifluidniNamesLts),
    ifluidnbNamesLts(_ifluidnbNamesLts), ifluidbNamesLts(_ifluidbNamesLts) 
{
  //---- stitch files from last time step
  
  // fluid
  stitchCGNS(fluidNamesLts, 0);
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
  transferStitchedToPartMb();

  std::cout << "RocRestartDriver created" << std::endl; 
} 

RocRestartDriver::~RocRestartDriver()
{
  for (int i = 0; i < stitchers.size(); ++i)
  {
    if (stitchers[i])
    { 
      delete stitchers[i]; // deletes mb and cg objs as well
      stitchers[i] = nullptr;
      mbObjs[i] = nullptr;
    }
  }
  for (int i = 0; i < fluidRmCg.size();++i)
  {
    if (fluidRmCg[i])
    {
      delete fluidRmCg[i];
      delete fluidRmMb[i];
      fluidRmCg[i] = nullptr;
      fluidRmMb[i] = nullptr;
    }
  } 
  for (int i = 0; i < ifluidNiRmCg.size();++i)
  {
    if (ifluidNiRmCg[i])
    {
      delete ifluidNiRmCg[i];
      delete ifluidNiRmMb[i];
      ifluidNiRmCg[i] = nullptr;
      ifluidNiRmMb[i] = nullptr;
    }
  } 
  for (int i = 0; i < ifluidBRmCg.size();++i)
  {
    if (ifluidBRmCg[i])
    {
      delete ifluidBRmCg[i];
      delete ifluidBRmMb[i];
      ifluidBRmCg[i] = nullptr;
      ifluidBRmMb[i] = nullptr;
    }
  } 
  for (int i = 0; i < ifluidNbRmCg.size();++i)
  {
    if (ifluidNbRmCg[i])
    {
      delete ifluidNbRmCg[i];
      delete ifluidNbRmMb[i];
      ifluidNbRmCg[i] = nullptr;
      ifluidNbRmMb[i] = nullptr;
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
  surfs.insert(surfs.begin(), mbObjs.begin()+1,mbObjs.end());
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

void RocRestartDriver::transferStitchedToPartMb()
{
  // fluid
  for (int i = 0; i < fluidRmMb.size(); ++i)
  {
    mbObjs[0]->transfer(fluidRmMb[i], "Consistent Interpolation"); 
  }
  // ifluid_ni
  for (int i = 0; i < ifluidNiRmMb.size(); ++i)
  {
    stitchedSurf->transfer(ifluidNiRmMb[i], "Consistent Interpolation");
  }
  // ifluid_nb
  std::cout << "nb size; " << ifluidNbRmMb.size() << std::endl;
  std::cout << ifluidnbNamesRm.size() << std::endl;
  for (int i = 0; i < ifluidNbRmMb.size(); ++i)
  {
    stitchedSurf->transfer(ifluidNbRmMb[i], "Consistent Interpolation");
  }
  // ifluid_b
  for (int i = 0; i < ifluidBRmMb.size(); ++i)
  {
    stitchedSurf->transfer(ifluidBRmMb[i], "Consistent Interpolation");
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
  std::string lastDir = inputjson["Rocout Last TS Directory"].as<std::string>();
  std::string base_t = inputjson["Rocout Last TS Base"].as<std::string>();

  std::vector<std::string> fluNamesRm(getCgFNames(prepDir, "fluid", "00.000000"));
  std::vector<std::string> ifluniNamesRm(getCgFNames(prepDir, "ifluid_ni", "00.000000"));
  std::vector<std::string> iflunbNamesRm(getCgFNames(prepDir, "ifluid_nb", "00.000000"));
  std::vector<std::string> iflubNamesRm(getCgFNames(prepDir, "ifluid_b", "00.000000"));
  std::vector<std::string> fluNamesLts(getCgFNames(lastDir, "fluid", base_t));
  std::vector<std::string> ifluniNamesLts(getCgFNames(lastDir, "ifluid_ni", base_t));
  std::vector<std::string> iflunbNamesLts(getCgFNames(lastDir, "ifluid_nb", base_t));
  std::vector<std::string> iflubNamesLts(getCgFNames(lastDir, "ifluid_b", base_t));
  
  RocRestartDriver* restartDriver 
    = new RocRestartDriver(fluNamesRm, ifluniNamesRm, iflunbNamesRm, iflubNamesRm,
                           fluNamesLts, ifluniNamesLts, iflunbNamesLts, iflubNamesLts);
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
