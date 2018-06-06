#include <RemeshDriver.H>
#include <meshStitcher.H>
#include <rocstarCgns.H>
#include <cgnsWriter.H>
#include <meshPartitioner.H>
#include <vtkAppendFilter.h>
#include <MeshGenDriver.H>
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
  #include <symmxParams.H>
#endif

// testing
#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGenericCell.h>
#include <vtkCell.h>
#include <vtkCellData.h>

RemeshDriver::RemeshDriver(const std::string& _case_dir, const std::string& _base_t,
                           const int _fluidproc, const int _ifluidniproc,
                           const int _ifluidbproc, const int _ifluidnbproc,
                           const json& remeshjson)
  : base_t(_base_t), case_dir(_case_dir), fluidproc(_fluidproc), ifluidniproc(_ifluidniproc),
    ifluidbproc(_ifluidbproc), ifluidnbproc(_ifluidnbproc)
{
  // stitch fluid files
  stitchFluidCGNS("fluid",fluidproc);
  // stitch ifluid ni files
  stitchIFluidCGNS(ifluidniName, "ifluid_ni", ifluidniproc);
  // stitch ifluid b files
  stitchIFluidCGNS(ifluidbName, "ifluid_b", ifluidbproc);
  // stitch ifluid nb files
  stitchIFluidCGNS(ifluidnbName, "ifluid_nb", ifluidnbproc);
  // get stitched meshes from stitcher vector
  for (int i = 0; i < stitchers.size(); ++i)
  {
    cgObjs.push_back(stitchers[i]->getStitchedCGNS());
    mbObjs.push_back(stitchers[i]->getStitchedMB());
  }
  // creates remeshedVol and remeshedSurf
  remesh(remeshjson);
  // creates stitchedSurf
  if (mbObjs.size() > 1)
  {
    stitchSurfaces();
    stitchedSurf->transfer(remeshedSurf, "Consistent Interpolation");
    stitchedSurf->write();
  }
  remeshedSurf->write();
  // partition the mesh
  partitionMesh();

  writeCobalt("tmp", "cobalt_test.cgr");
  std::cout << "RemeshDriver created" << std::endl;
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
  if (mshgendrvr)
  {
    delete mshgendrvr;
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
  std::cout << "RemeshDriver destroyed" << std::endl;
}

void RemeshDriver::stitchFluidCGNS(const std::string& prefix, const int numproc)
{
  if (numproc)
  {
    setCgFnames(fluidNames, prefix, numproc);
    stitchers.push_back(new meshStitcher(fluidNames));
  }
}

void RemeshDriver::stitchIFluidCGNS(std::string& name, const std::string& prefix, const int numproc)
{
  if (numproc)
  {
    setCgFname(name, prefix);
    stitchers.push_back(new meshStitcher(numproc, name));
  }
} 

void RemeshDriver::remesh(const json& remeshjson)
{
  std::cout << "Extracting surface mesh #############################################\n";
  std::unique_ptr<meshBase> surf 
    = std::unique_ptr<meshBase>(meshBase::Create(mbObjs[0]->extractSurface(), "skinMeshPart.vtp")); 
  // remeshing with engine specified in input
  surf->write("skinMeshPart.stl"); 
  mshgendrvr = MeshGenDriver::readJSON("skinMeshPart.stl", "remeshedVol.vtu", remeshjson);
  remeshedVol = mshgendrvr->getNewMesh();
  remeshedSurf = meshBase::Create(remeshedVol->extractSurface(), "remeshedSurf.vtp");
}

void RemeshDriver::stitchSurfaces()
{
  // stitch b, ni and nb surfaces
  vtkSmartPointer<vtkAppendFilter> appender
    = vtkSmartPointer<vtkAppendFilter>::New();
  for (int i = 1; i < mbObjs.size(); ++i)
  {
    appender->AddInputData(mbObjs[i]->getDataSet());
  }
  appender->Update();
  stitchedSurf = meshBase::Create(appender->GetOutput(), "stitchedSurf.vtu");
  stitchedSurf->setContBool(0);
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
  RemeshDriver* remeshdrvobj = new RemeshDriver(case_dir, base_t, fluidproc, ifluidniproc, 
                                                ifluidbproc, ifluidnbproc, remeshjson);
  return remeshdrvobj;
}

std::map<int,int> RemeshDriver::readPatchMap(const std::string& mapFile)
{
  std::ifstream inputStream(mapFile); 
  if (!inputStream.good())
  {
    std::cout << "error opening file " << mapFile << std::endl;
    exit(1);
  }
  
  int numPatch;
  int numMap;
  std::string line;
  getline(inputStream,line);
  std::stringstream ss(line);
  ss >> numPatch;
  getline(inputStream,line);
  std::stringstream ss1(line);
  ss1 >> numMap;
  std::map<int, int> patchMap; 
  while(getline(inputStream, line))
  {
    std::stringstream ss(line);
    int lo; 
    int up; 
    int patch;
    ss >> lo >> up >> patch;
    for (int i = lo; i <= up; ++i)
    {   
      patchMap[i] = patch;
    }   
  }   
  return patchMap;
}


void RemeshDriver::writeCobalt(const std::string& mapFile, std::ofstream& outputStream)
{

  if (!remeshedSurf) 
  {
    std::cout << "remeshed surface has not been generated!" << std::endl;
    exit(1);
  }
  
  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New(); 
  vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New(); 
  std::map<std::vector<int>, std::pair<int,int>> faceMap;

  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkCellLocator> surfCellLocator = remeshedSurf->buildLocator(); 
  // maximum number of vertices per face (to be found in proceeding loop)
  int nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0; 

  for (int i = 0; i < remeshedVol->getNumberOfCells(); ++i)
  {
    // get cell i 
    remeshedVol->getDataSet()->GetCell(i, genCell1);
    // get faces, find cells sharing it. if no cell shares it, 
    // use the locator of the remeshedSurf to find the patch number
    int numFaces = genCell1->GetNumberOfFaces();
    nFacesPerCellMax = (nFacesPerCellMax < numFaces ? numFaces : nFacesPerCellMax);
    for (int j = 0; j < numFaces; ++j)
    {
      vtkCell* face = genCell1->GetFace(j);
      bool shared = 0;
      int numVerts = face->GetNumberOfPoints();
      nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
      facePtIds = face->GetPointIds(); 
      remeshedVol->getDataSet()->GetCellNeighbors(i, facePtIds, sharedCellPtIds); 
      std::vector<int> facePntIds(numVerts);
      for (int k = 0; k < numVerts; ++k)
      {
        facePntIds[k] = face->GetPointId(k)+1;
      }
      if (sharedCellPtIds->GetNumberOfIds())
      {

        faceMap.insert(std::pair<std::vector<int>, std::pair<int,int>> 
                        (facePntIds, std::make_pair(i+1, (int) sharedCellPtIds->GetId(0))));
      }
      else
      {
        int id = face->GetPointId(0);
        double x[3];
        vtkIdType closestCellId;
        int subId;
        double minDist2;
        double closestPoint[3];
        remeshedVol->getDataSet()->GetPoint(id,x);
        // find closest point and closest cell to x
        surfCellLocator->FindClosestPoint(x, closestPoint, genCell2,closestCellId,subId,minDist2);
        double patchNo[1];
        remeshedSurf->getDataSet()->GetCellData()->GetArray("patchNo")
                                                  ->GetTuple(closestCellId, patchNo);

        faceMap.insert(std::pair<std::vector<int>, std::pair<int,int>>
                        (facePntIds, std::make_pair(i+1, (int) -1*patchNo[0])));

      }
    }
  }
  
  std::map<int,int> patchMap;
  for (int i = 0; i < remeshedSurf->getNumberOfCells(); ++i)
  {
    double patchNo[1];
    remeshedSurf->getDataSet()->GetCellData()->GetArray("patchNo")
                                              ->GetTuple(i, patchNo);
    patchMap.insert(std::pair<int,int>(patchNo[0],i));
  }
   
  // write cobalt file
  outputStream << 3 << " " << 1 << " " << patchMap.size() << std::endl;
  outputStream << remeshedVol->getNumberOfPoints() << " " << faceMap.size()
               << " " << remeshedVol->getNumberOfCells() << " "
               << nVerticesPerFaceMax << " " << nFacesPerCellMax << std::endl;
  for (int i = 0; i < remeshedVol->getNumberOfPoints(); ++i)
  {
    std::vector<double> pnt(remeshedVol->getPoint(i));
    outputStream << "  "  << std::setw(21) << std::setprecision(15)
                 << pnt[0] << "   " << pnt[1] << "   " << pnt[2] << std::endl;
  }

  auto it = faceMap.begin();
  while (it != faceMap.end())
  {
    outputStream << it->first.size() << " ";
    for (int i = 0; i < it->first.size(); ++i)
    {
      outputStream << it->first[i] << " ";
    }
    outputStream << it->second.first << " " << it->second.second << std::endl;
    ++it;
  }
}

void RemeshDriver::writeCobalt(const std::string& mapFile, const std::string& ofname)
{
  std::ofstream outputStream(ofname);
  if(!outputStream.good()) 
  {
    std::cout << "Cannot open file " << ofname << std::endl;
    exit(1);
  }
  writeCobalt(mapFile,outputStream); 
}
