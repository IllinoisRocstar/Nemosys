#include <RemeshDriver.H>
#include <meshStitcher.H>
#include <rocstarCgns.H>
#include <cgnsWriter.H>
#include <meshPartitioner.H>
#include <vtkAppendFilter.h>
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
  #include <symmxParams.H>
#endif

RemeshDriver::RemeshDriver(const std::string& _case_dir, const std::string& _base_t,
                           const int _fluidproc, const int _ifluidniproc,
                           const int _ifluidbproc, const int _ifluidnbproc,
                           meshingParams* _params)
  : base_t(_base_t), case_dir(_case_dir), params(_params), fluidproc(_fluidproc),
    ifluidniproc(_ifluidniproc), ifluidbproc(_ifluidbproc), ifluidnbproc(_ifluidnbproc)
{
  // stitch fluid files
  if (fluidproc)
  {
    setCgFnames(fluidNames, "fluid", fluidproc);
    stitchers.push_back(new meshStitcher(fluidNames));
  }
  // stitch ifluid ni files
  if (ifluidniproc)
  {
    setCgFname(ifluidniName, "ifluid_ni");
    stitchers.push_back(new meshStitcher(ifluidniproc, ifluidniName));
  }  
  // stitch ifluid b files
  if (ifluidbproc)
  {
    setCgFname(ifluidbName, "ifluid_b");
    stitchers.push_back(new meshStitcher(ifluidbproc, ifluidbName));
  }
  // stitch ifluid nb files
  if (ifluidnbproc)
  {
    setCgFname(ifluidnbName, "ifluid_nb");
    stitchers.push_back(new meshStitcher(ifluidnbproc, ifluidnbName));
  }

  for (int i = 0; i < stitchers.size(); ++i)
  {
    cgObjs.push_back(stitchers[i]->getStitchedCGNS());
    mbObjs.push_back(stitchers[i]->getStitchedMB());
  }

  // remesh the volume with symmetrix
  symmxRemesh();
  
  // stitch b, ni and nb surfaces
  vtkSmartPointer<vtkAppendFilter> appender
    = vtkSmartPointer<vtkAppendFilter>::New();
  appender->AddInputData(mbObjs[1]->getDataSet());
  appender->AddInputData(mbObjs[2]->getDataSet());
  appender->Update();
  stitchedSurf = meshBase::Create(appender->GetOutput(), "stitchedSurf.vtu");
  stitchedSurf->transfer(remeshedSurf, "Consistent Interpolation");
  stitchedSurf->write();
  remeshedSurf->write();
  // partition the mesh
  partitionMesh();
}

RemeshDriver::~RemeshDriver()
{
  for (int i = 0; i < stitchers.size(); ++i)
  {
    if (stitchers[i])
    { 
      delete stitchers[i]; // deletes mb and cg objs as well
      stitchers[i] = nullptr;
      cgObjs[i] = nullptr;
      mbObjs[i] = nullptr;
    }
  }
  if (remeshedVol)
  {
    delete remeshedVol;
    remeshedVol = nullptr;
  }
  if (remeshedSurf)
  {
    delete remeshedSurf;
    remeshedSurf = nullptr;
  }
  if (stitchedSurf)
  {
    delete stitchedSurf;
    stitchedSurf = nullptr;
  } 
  if (params)
  {
    delete params;
    params = nullptr;
  }
}

void RemeshDriver::symmxRemesh()
{
  std::cout << "Extracting surface mesh #############################################\n";
  meshBase* surf = meshBase::Create(mbObjs[0]->extractSurface(), "skinMeshPart.vtp"); 
  // remeshing with simmetrix
  surf->write("skinMeshPart.stl"); 
  remeshedVol = meshBase::generateMesh("skinMeshPart.stl","simmetrix",params);
  mbObjs[0]->setContBool(0);
  remeshedSurf = meshBase::Create(remeshedVol->extractSurface(), "remeshedSurf.vtp");
  delete surf; surf = nullptr;
}

void RemeshDriver::partitionMesh()
{
  std::cout << " Partitioning the mesh with METIS ###################################\n";
  meshPartitioner* mPart = new meshPartitioner(remeshedVol);
  mPart->partition(fluidproc);
  
  // write CGNS files for the new grid
  for (int iCg=0; iCg < fluidproc; ++iCg)
  {
    std::string fCgName(fluidNames[iCg]);
    std::size_t pos = fCgName.find_last_of("/");
    fCgName = fCgName.substr(pos+1);
    fCgName = trim_fname(fCgName,"New.cgns");
    std::cout << "Writing remeshed cgns part to " << fCgName << std::endl;
    // define elementary information
    cgnsWriter* cgWrtObj = new cgnsWriter(fCgName, cgObjs[0]->getBaseName(), 3, 3);
    cgWrtObj->setUnits(cgObjs[0]->getMassUnit(), cgObjs[0]->getLengthUnit(),
		 cgObjs[0]->getTimeUnit(), cgObjs[0]->getTemperatureUnit(),
		 cgObjs[0]->getAngleUnit());
    cgWrtObj->setBaseItrData(cgObjs[0]->getBaseItrName(), 
                             cgObjs[0]->getNTStep(), 
                             cgObjs[0]->getTimeStep());
    cgWrtObj->setZoneItrData(cgObjs[0]->getZoneItrName(), 
                             cgObjs[0]->getGridCrdPntr(), 
                             cgObjs[0]->getSolutionPntr());
    cgWrtObj->setZone(cgObjs[0]->getZoneName(iCg), cgObjs[0]->getZoneType());
    cgWrtObj->setNVrtx(mPart->getNNdePart(iCg));
    cgWrtObj->setNCell(mPart->getNElmPart(iCg));
    // define coordinates
    std::vector<std::vector<double>> comp_crds(remeshedVol->getVertCrds()); 
    cgWrtObj->setGridXYZ(mPart->getCrds(iCg, comp_crds[0]), 
		                     mPart->getCrds(iCg, comp_crds[1]), 
		                     mPart->getCrds(iCg, comp_crds[2]));
    // define connctivity
    cgWrtObj->setSection(cgObjs[0]->getSectionName(), 
		                     (ElementType_t) cgObjs[0]->getElementType(), 
		                     mPart->getConns(iCg));

    // write partitioned vtk files with data transfered from stitched mesh
    std::vector<int> vtkConn(mPart->getConns(iCg));
    for (auto it = vtkConn.begin(); it != vtkConn.end(); ++it)
    {
      *it -= 1;
    }
    std::stringstream vtkname;
    vtkname << "partition" << iCg << ".vtu";
    meshBase* mbPart = meshBase::Create(mPart->getCrds(iCg, comp_crds[0]),
                                        mPart->getCrds(iCg, comp_crds[1]),
                                        mPart->getCrds(iCg, comp_crds[2]),
                                        vtkConn, VTK_TETRA, vtkname.str());
    mbObjs[0]->transfer(mbPart, "Consistent Interpolation");
    mbPart->write();

    // define vertex and cell data 
    std::map<std::string, GridLocation_t> slnNLMap = cgObjs[0]->getSolutionNameLocMap();
    for (auto is=slnNLMap.begin(); is!=slnNLMap.end(); is++)
      cgWrtObj->setSolutionNode(is->first, is->second);
    // write skeleton of the file
    cgWrtObj->writeGridToFile();

    /////////////////////// WRITE SOLUTION DATA TO CGNS ///////////////////////////////

    // write individual data fields
    std::map<int,std::pair<int,keyValueList> > slnMap = cgObjs[0]->getSolutionMap();
    std::vector<GridLocation_t> gLoc = cgObjs[0]->getSolutionGridLocations();
    std::vector<std::string> slnName = cgObjs[0]->getSolutionNodeNames();

    int iSol = -1;
    for (auto is=slnMap.begin(); is!=slnMap.end(); is++)
    {
      std::pair<int,keyValueList> slnPair = is->second;
      int slnIdx = slnPair.first;
      keyValueList fldLst = slnPair.second;
      for (auto ifl=fldLst.begin(); ifl!=fldLst.end(); ifl++)
      {
	      iSol++;
	      std::vector<double> partPhysData;
	      int nData;
	      if (gLoc[iSol] == Vertex)
	      {
	       mbPart->getPointDataArray(ifl->second, partPhysData);              
	      } 
        else 
        {
          mbPart->getCellDataArray(ifl->second, partPhysData);
	      }
	      std::cout << "Writing "
	          << partPhysData.size() 
	          << " to "
	          << ifl->second
	          << " located in "
	          << slnName[iSol]
	          << std::endl;
	      // write to file
	      cgWrtObj->writeSolutionField(ifl->second, slnName[iSol], RealDouble, &partPhysData[0]);
      }
    }
    delete cgWrtObj; cgWrtObj = nullptr;
    delete mbPart; mbPart = nullptr;
  }
  delete mPart; mPart = nullptr;
}


void RemeshDriver::setCgFnames(std::vector<std::string>& names, const std::string& prefix,
                               const int numproc)
{
  names.resize(numproc);
  for (int i = 0; i < numproc; ++i)
  {
    std::stringstream basename;
    basename << case_dir << "/" << prefix << "_" << base_t << "_";
    if (i < 10)
      basename << "000" << i << ".cgns";
    else if (i >= 10 && i < 100)
      basename << "00" << i << ".cgns";
    else if (i >= 100 && i < 1000)
      basename << 0 << i << ".cgns";
    names[i] = basename.str();
  } 
}

void RemeshDriver::setCgFname(std::string& name, const std::string& prefix)
{
  std::stringstream basename;
  basename << case_dir << "/" << prefix << "_" << base_t << "_0000.cgns";
  name = basename.str();
}

RemeshDriver* RemeshDriver::readJSON(json inputjson)
{
  std::string case_dir = inputjson["Case Directory"].as<std::string>();
  std::string base_t = inputjson["Base Time Step"].as<std::string>();
  json numprocjson = inputjson["Number Of Processors"];
  int fluidproc = numprocjson["fluid"].as<int>();
  int ifluidniproc = numprocjson["ifluid_ni"].as<int>();
  int ifluidbproc = numprocjson["ifluid_b"].as<int>();
  int ifluidnbproc = numprocjson["ifluid_nb"].as<int>();

  json remeshjson = inputjson["Remeshing Options"];

  std::string meshEngine = remeshjson["Mesh Generation Engine"].as<std::string>();

  if (!meshEngine.compare("simmetrix"))
  {
    #ifndef HAVE_SYMMX
      std::cerr << "Nemosys must be recompiled with simmetrix support" << std::endl;
      exit(1);
    #else
      if (!remeshjson.has_key("License File"))
      {
        std::cerr << "Simmetrix License file must be specified in json" << std::endl;
        exit(1);
      }
      
      symmxParams* params = new symmxParams();
      params->licFName = remeshjson["License File"].as<std::string>();
      params->features = remeshjson["Features"].as<std::string>();
      params->logFName = remeshjson["Log File"].as<std::string>();
      
      std::string defaults = remeshjson["Meshing Parameters"]["Simmetrix Parameters"].as<std::string>();
      if (!defaults.compare("default"))
      {
        RemeshDriver* remeshdrvobj = new RemeshDriver(case_dir, base_t, fluidproc,
                                                      ifluidniproc, ifluidbproc, ifluidnbproc,  
                                                      params);     
        return remeshdrvobj;
      }
      else
      {
        json symmxparams = remeshjson["Meshing Parameters"]["Simmetrix Parameters"];
        if (symmxparams.has_key("Mesh Size"))
          params->meshSize = symmxparams["Mesh Size"].as<double>();
        if (symmxparams.has_key("Anisotropic Curvature Refinement"))
          params->anisoMeshCurv = symmxparams["Anisotropic Curvature Refinement"].as<double>();
        if (symmxparams.has_key("Global Gradation Rate"))
          params->glbSizeGradRate = symmxparams["Global Gradation Rate"].as<double>();
        if (symmxparams.has_key("Surface Mesh Improver Gradation Rate"))
          params->surfMshImprovGradRate = symmxparams["Surface Mesh Improver Gradation Rate"].as<double>();
        if (symmxparams.has_key("Surface Mesh Improver Min Size"))
          params->surfMshImprovMinSize = symmxparams["Surface Mesh Improver Min Size"].as<double>();
      
        RemeshDriver* remeshdrvobj = new RemeshDriver(case_dir, base_t, fluidproc,
                                                      ifluidniproc, ifluidbproc, ifluidnbproc,  
                                                      params);     
        return remeshdrvobj;
      } 
    #endif
  }

  else
  {
    std::cout << "Mesh generation engine " << meshEngine << " is not supported" << std::endl;
    exit(1);
  }
}
