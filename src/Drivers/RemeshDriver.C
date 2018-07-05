// nemosys
#include <RemeshDriver.H>
#include <meshStitcher.H>
#include <rocstarCgns.H>
#include <cgnsWriter.H> // not used
#include <meshPartitioner.H> // not used
#include <MeshGenDriver.H>
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
                           const json& remeshjson)
  : fluidNames(_fluidNames), ifluidniNames(_ifluidniNames), ifluidnbNames(_ifluidnbNames),
    ifluidbNames(_ifluidbNames)
{
  // stitch fluid files
  stitchCGNS(fluidNames,0);
  // stitch ifluid_ni files
  stitchCGNS(ifluidniNames,1);
  // stitch ifluid_b files
  stitchCGNS(ifluidbNames,1);
  // stitch ifluid_ng files
  stitchCGNS(ifluidnbNames,1);
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
    std::vector<std::string> patchno(1,"patchNo");
    stitchedSurf->transfer(remeshedSurf, "Consistent Interpolation", patchno, 1);
    stitchedSurf->write();
  }
  remeshedSurf->write();
  writeCobalt("remeshedVol.cgi", "remeshedVol.cgr");
  std::vector<meshBase*> mbPartitions(meshBase::partition(remeshedVol, 10));
  for (int i = 0; i < mbPartitions.size(); ++i)
  {
    mbPartitions[i]->write();
    delete mbPartitions[i];
  }
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
    mshgendrvr = nullptr;
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


void RemeshDriver::stitchCGNS(const std::vector<std::string>& fnames, bool surf)
{
  if (fnames.size())
  {
    stitchers.push_back(new meshStitcher(fnames, surf));
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
  mshgendrvr = MeshGenDriver::readJSON("extractedSurface.stl", "remeshedVol.vtu", remeshjson);
  remeshedVol = mshgendrvr->getNewMesh();
  remeshedSurf = meshBase::Create(remeshedVol->extractSurface(), "remeshedSurf.vtp");
}

void RemeshDriver::stitchSurfaces()
{
  // stitch b, ni and nb surfaces
  std::vector<meshBase*> surfs;
  surfs.insert(surfs.begin(), mbObjs.begin()+1,mbObjs.end());
  stitchedSurf = meshBase::stitchMB(surfs);
  stitchedSurf->setContBool(0);
  stitchedSurf->setFileName("stitchedSurf.vtu");
}

RemeshDriver* RemeshDriver::readJSON(json inputjson)
{
  std::string case_dir = inputjson["Case Directory"].as<std::string>();
  std::string base_t = inputjson["Base Time Step"].as<std::string>();
  std::vector<std::string> fluNames(getCgFNames(case_dir, "fluid", base_t));
  std::vector<std::string> ifluniNames(getCgFNames(case_dir, "ifluid_ni", base_t));
  std::vector<std::string> iflunbNames(getCgFNames(case_dir, "ifluid_nb", base_t));
  std::vector<std::string> iflubNames(getCgFNames(case_dir, "ifluid_b", base_t));
  json remeshjson = inputjson["Remeshing Options"];
  RemeshDriver* remeshdrvobj = new RemeshDriver(fluNames, ifluniNames, 
                                                iflunbNames, iflubNames, remeshjson);
  return remeshdrvobj;
}

void RemeshDriver::writePatchMap(const std::string& mapFile, const std::map<int,int>& patchMap)
{
  std::ofstream outputStream(mapFile);
  if (!outputStream.good())
  {
    std::cout << "Error opening file " << mapFile << std::endl;
    exit(1);
  }
  writePatchMap(outputStream, patchMap);
}

void RemeshDriver::writePatchMap(std::ofstream& outputStream, const std::map<int,int>& patchMap)
{
  outputStream << patchMap.size() << std::endl;
  outputStream << patchMap.size() << std::endl;
  auto it = patchMap.begin();
  while (it != patchMap.end())
  {
    for (int i = 0; i < 3; ++i)
    {
      outputStream << std::setw(2) << std::left << it->first << " ";
    }
    outputStream << std::endl;
    ++it;
  }
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
  std::map<std::vector<int>, std::pair<int,int>, sortIntVec_compare> faceMap;
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
                        (facePntIds, std::make_pair(i+1, (int) sharedCellPtIds->GetId(0)+1))); 
      }
      else
      {
        double p1[3], p2[3], p3[3];
        face->GetPoints()->GetPoint(0,p1);
        face->GetPoints()->GetPoint(1,p2);
        face->GetPoints()->GetPoint(2,p3);
        double faceCenter[3];
        for (int k = 0; k < numVerts; ++k)
        {
          faceCenter[k] = (p1[k] + p2[k] + p3[k])/3.0;
        } 
        vtkIdType closestCellId;
        int subId;
        double minDist2;
        double closestPoint[3];
        // find closest point and closest cell to faceCenter
        surfCellLocator->FindClosestPoint(faceCenter, closestPoint, genCell2,closestCellId,subId,minDist2);
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

  // write patch mapping file
  writePatchMap(mapFile, patchMap);
  // write cobalt file
  outputStream << 3 << "   " << 1 << "  " << patchMap.size() << std::endl;
  outputStream << remeshedVol->getNumberOfPoints() << " " << faceMap.size()
               << " " << remeshedVol->getNumberOfCells() << " "
               << nVerticesPerFaceMax << " " << nFacesPerCellMax << std::endl;
  for (int i = 0; i < remeshedVol->getNumberOfPoints(); ++i)
  {
    std::vector<double> pnt(remeshedVol->getPoint(i));
    outputStream << std::setw(21) << std::fixed << std::setprecision(15)
                 << pnt[0] << "   " << pnt[1] << "   " << pnt[2] << std::endl;
  }
  
  auto it = faceMap.begin();
  while (it != faceMap.end())
  {
    outputStream << it->first.size() << " ";
    auto faceIdIter = it->first.begin();
    while (faceIdIter != it->first.end())
    {
      outputStream << *faceIdIter << " ";
      ++faceIdIter;
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

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
}

bool sortIntVec_compare::operator() (std::vector<int> lhs, 
                                     std::vector<int> rhs) const
{
  std::sort(lhs.begin(), lhs.end());
  std::sort(rhs.begin(), rhs.end());
  return lhs < rhs;
}

//void RemeshDriver::partitionMesh()
//{
//  std::cout << " Partitioning the mesh with METIS ###################################\n";
//  meshPartitioner* mPart = new meshPartitioner(remeshedVol);
////  mPart->partition(fluidNames.size());
//  mPart->partition(10);
//  // write CGNS files for the new grid
//  for (int iCg=0; iCg < fluidNames.size(); ++iCg)
//  {
//    std::size_t pos = fluidNames[iCg].find_last_of("/");
//    std::string fCgName(fluidNames[iCg].substr(pos+1));
//    fCgName = trim_fname(fCgName,"New.cgns");
//    std::cout << "Writing remeshed cgns part to " << fCgName << std::endl;
//    // define elementary information
//    cgnsWriter* cgWrtObj = new cgnsWriter(fCgName, cgObjs[0]->getBaseName(), 3, 3);
//    cgWrtObj->setUnits(cgObjs[0]->getMassUnit(), cgObjs[0]->getLengthUnit(),
//		 cgObjs[0]->getTimeUnit(), cgObjs[0]->getTemperatureUnit(),
//		 cgObjs[0]->getAngleUnit());
//    cgWrtObj->setBaseItrData(cgObjs[0]->getBaseItrName(), 
//                             cgObjs[0]->getNTStep(), 
//                             cgObjs[0]->getTimeStep());
//    cgWrtObj->setZoneItrData(cgObjs[0]->getZoneItrName(), 
//                             cgObjs[0]->getGridCrdPntr(), 
//                             cgObjs[0]->getSolutionPntr());
//    cgWrtObj->setZone(cgObjs[0]->getZoneName(iCg), cgObjs[0]->getZoneType());
//    cgWrtObj->setNVrtx(mPart->getNNdePart(iCg));
//    cgWrtObj->setNCell(mPart->getNElmPart(iCg));
//    // define coordinates
//    std::vector<std::vector<double>> comp_crds(remeshedVol->getVertCrds()); 
//    cgWrtObj->setGridXYZ(mPart->getCrds(iCg, comp_crds[0]), 
//		                     mPart->getCrds(iCg, comp_crds[1]), 
//		                     mPart->getCrds(iCg, comp_crds[2]));
//    // define connctivity
//    cgWrtObj->setSection(cgObjs[0]->getSectionName(), 
//		                     (ElementType_t) cgObjs[0]->getElementType(), 
//		                     mPart->getConns(iCg));
//
//    // write partitioned vtk files with data transfered from stitched mesh
//    std::vector<int> vtkConn(mPart->getConns(iCg));
//    for (auto it = vtkConn.begin(); it != vtkConn.end(); ++it)
//    {
//      *it -= 1;
//    }
//    std::stringstream vtkname;
//    vtkname << "partition" << iCg << ".vtu";
//    meshBase* mbPart = meshBase::Create(mPart->getCrds(iCg, comp_crds[0]),
//                                        mPart->getCrds(iCg, comp_crds[1]),
//                                        mPart->getCrds(iCg, comp_crds[2]),
//                                        vtkConn, VTK_TETRA, vtkname.str());
//    mbObjs[0]->transfer(mbPart, "Consistent Interpolation");
//    mbPart->write();
//
//    // define vertex and cell data 
//    std::map<std::string, GridLocation_t> slnNLMap = cgObjs[0]->getSolutionNameLocMap();
//    for (auto is=slnNLMap.begin(); is!=slnNLMap.end(); is++)
//      cgWrtObj->setSolutionNode(is->first, is->second);
//    // write skeleton of the file
//    cgWrtObj->writeGridToFile();
//
//    /////////////////////// WRITE SOLUTION DATA TO CGNS ///////////////////////////////
//
//    // write individual data fields
//    std::map<int,std::pair<int,keyValueList> > slnMap = cgObjs[0]->getSolutionMap();
//    std::vector<GridLocation_t> gLoc = cgObjs[0]->getSolutionGridLocations();
//    std::vector<std::string> slnName = cgObjs[0]->getSolutionNodeNames();
//
//    int iSol = -1;
//    for (auto is=slnMap.begin(); is!=slnMap.end(); is++)
//    {
//      std::pair<int,keyValueList> slnPair = is->second;
//      int slnIdx = slnPair.first;
//      keyValueList fldLst = slnPair.second;
//      for (auto ifl=fldLst.begin(); ifl!=fldLst.end(); ifl++)
//      {
//	      iSol++;
//	      std::vector<double> partPhysData;
//	      int nData;
//	      if (gLoc[iSol] == Vertex)
//	      {
//	       mbPart->getPointDataArray(ifl->second, partPhysData);              
//	      } 
//        else 
//        {
//          mbPart->getCellDataArray(ifl->second, partPhysData);
//	      }
//	      std::cout << "Writing "
//	          << partPhysData.size() 
//	          << " to "
//	          << ifl->second
//	          << " located in "
//	          << slnName[iSol]
//	          << std::endl;
//	      // write to file
//	      cgWrtObj->writeSolutionField(ifl->second, slnName[iSol], RealDouble, &partPhysData[0]);
//      }
//    }
//    delete cgWrtObj; cgWrtObj = nullptr;
//    delete mbPart; mbPart = nullptr;
//  }
//  delete mPart; mPart = nullptr;
//}
