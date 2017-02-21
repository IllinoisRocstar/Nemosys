/* Special purpose class for Rocstar CGNS files */


#include <gridTransfer.H>

///////////////////////////////////////////////////
// INITIALIZATION
///////////////////////////////////////////////////

gridTransfer::gridTransfer(std::string srcFname, std::string trgFname) :
              cgnsAnalyzer(srcFname), 
              srcCgFName(srcFname), trgCgFName(trgFname),
              srcModel(NULL), trgModel(NULL),
              srcMesh(NULL), trgMesh(NULL)
{
  // source CGNS file processing
  std::size_t _loc = srcCgFName.find_last_of("_");
  baseCgFNameSrc = srcCgFName.substr(0,_loc+1);
  padSizeSrc = srcCgFName.size() - baseCgFNameSrc.size() - 5;

  // target CGNS file processing
  _loc = trgCgFName.find_last_of("_");
  baseCgFNameTrg = trgCgFName.substr(0,_loc+1);
  padSizeTrg = trgCgFName.size() - baseCgFNameTrg.size() - 5;
}

gridTransfer::~gridTransfer()
{
 // cleaning up
 for (int ic=0; ic<srcCgObjs.size(); ic++)
   delete srcCgObjs[ic];
 for (int ic=0; ic<trgCgObjs.size(); ic++)
   delete trgCgObjs[ic];
}

/*
   Loads a series of source cgns files based on current file name
*/
void gridTransfer::loadSrcCgSeries(int nCg)
{
  for (int iCg=0; iCg<nCg; iCg++)
  {
    std::ostringstream suffix1;
    std::ostringstream suffix2;
    suffix1 << iCg;
    std::string tmp = suffix1.str();
    for (int iPad=0; iPad<padSizeSrc-tmp.size(); iPad++)
      suffix2 << "0";
    suffix2 << tmp << ".cgns";
    std::string fName = baseCgFNameSrc + suffix2.str();
    cgnsAnalyzer* cgTmp = new cgnsAnalyzer(fName);
    // load base grid
    cgTmp->loadGrid();
    srcCgObjs.push_back(cgTmp);
    srcCgFNames.push_back(fName);
    // loading all zones
    for (int iz = 1; iz<=srcCgObjs[iCg]->getNZone(); iz++)
    {
      std::cout << "  Zone ID : " << getZoneName(srcCgObjs[iCg], iz) << std::endl;
      stitchMe(srcCgObjs[iCg], iz);
      std::cout << "Finished processing of the zone " << getZoneName(srcCgObjs[iCg], iz) << std::endl;
    }
  }

}

/*
   Loads a series of target cgns files based on current file name
*/
void gridTransfer::loadTrgCg()
{
  cgnsAnalyzer* cgTmp = new cgnsAnalyzer(trgCgFName);
  // loading base grid
  cgTmp->loadGrid();
  trgCgObjs.push_back(cgTmp);
  trgCgFNames.push_back(trgCgFName);
}

////////////////////////////////////////////////////
// DATA PROCESSING
////////////////////////////////////////////////////
void gridTransfer::gridStats()
{
  std::cout << "  Number of Points    : " << srcMesh->nbPoints << std::endl;
  std::cout << "  Number of Edges     : " << srcMesh->nbEdges << std::endl;
  std::cout << "  Number of Triangles : " << srcMesh->nbTriangles << std::endl;
  std::cout << "  Number of Quads     : " << srcMesh->nbQuads << std::endl;
  std::cout << "  Number of Tets      : " << srcMesh->nbTets << std::endl;
  std::cout << "  Number of Hexes     : " << srcMesh->nbHexes << std::endl;
  std::cout << "  Number of Prisms    : " << srcMesh->nbPrisms << std::endl;

  // running checkMesh
  //MAd::checkMesh(srcMesh);

  // running adapter checks on the grid
  classifyMAdMeshBnd(srcMesh);
  MAd::PWLSField* sizeField = new MAd::PWLSField(srcMesh);
  sizeField->setCurrentSize();
  MAd::MeshAdapter* ma;
  ma = new MAd::MeshAdapter(srcMesh, sizeField);
  ma->printStatistics(std::cout);
  delete(ma);
}

void gridTransfer::gridCheck()
{
  // running checkMesh
  MAd::checkMesh(srcMesh);
}

void gridTransfer::transfer()
{
  // exporting to MAd if not yet
  if (!srcMesh)
    exportMeshToMAdLib("src");
  if (!trgMesh)
    exportMeshToMAdLib("trg");

  // transfering solution values to the new mesh
  std::vector<std::string> slnList;
  getSolutionDataNames(slnList);
  solution_type_t st;
  std::vector<double> srcVrtCrds = getVertexCoords();
  std::vector<double> srcElmCntCrds = getElmCntCoords(srcMesh);
  std::vector<double> trgElmCntCrds = getElmCntCoords(trgMesh);

  // preparing interpolators
  basicInterpolant interpNde = basicInterpolant(3, nVertex, 4, srcVrtCrds);
  basicInterpolant interpElm = basicInterpolant(3, nElem, 4, srcElmCntCrds);
  // loop on solutions
  for (auto is=slnList.begin(); is!=slnList.end(); is++)
  {
    std::vector<double> srcSlnVec;
    std::vector<double> trgSlnVec;
    int badPnt=0;
    int outNData, outNDim;
    int inNData = 0;
    st = getSolutionDataStitched(*is, srcSlnVec, outNData, outNDim);
    if (st == NODAL)
    {
      // nodal value transfer
      std::cout << "Transfering nodal " << *is << std::endl;
      for (int iNde=0; iNde<trgCgObjs[0]->getNVertex(); iNde++)
      {
	std::vector<double> vrtCrds, prms;
	std::vector<int>  vrtIds;
	vrtCrds = trgCgObjs[0]->getVertexCoords(iNde);
	//std::cout << vrtCrds[0] << " " << vrtCrds[1] << " " << vrtCrds[2] << std::endl;
	int elmIdx = getBaryCrds(vrtCrds, prms, vrtIds);
	if (elmIdx < 0)
	{
          badPnt++;
	  //std::cout << badPnt << std::endl;
	  std::vector<double> trgVrtData;
	  interpNde.interpolate(1, vrtCrds, srcSlnVec, trgVrtData);
	  trgSlnVec.push_back(trgVrtData[0]);
	}
	else
	  trgSlnVec.push_back( prms[0]*srcSlnVec[vrtIds[0]] +
			       prms[1]*srcSlnVec[vrtIds[1]] +
			       prms[2]*srcSlnVec[vrtIds[2]] +
			       prms[3]*srcSlnVec[vrtIds[3]] );
      }
      inNData = trgCgObjs[0]->getNVertex();
      std::cout << "Finished transfering with " << badPnt << " bad nodes." << std::endl;
    }
    else if (st == ELEMENTAL)
    {
      // elemental value transfer
      std::cout << "Transfering elemental " << *is << std::endl;
      interpElm.interpolate(trgCgObjs[0]->getNElement(), trgElmCntCrds, srcSlnVec, trgSlnVec);
      inNData = trgCgObjs[0]->getNElement();
    }
    trgCgObjs[0]->appendSolutionData(*is, trgSlnVec, st, inNData , 3);
  }
}

void gridTransfer::writeTrgCg(std::string cgFName)
{
  std::cout << "Writing remeshed " << cgFName << std::endl;
  cgnsAnalyzer* cgObj1 = trgCgObjs[0];
  // define elementary information
  cgnsWriter* cgWrtObj = new cgnsWriter(cgFName, cgObj1->getBaseName(), 3, 3);
  cgWrtObj->setUnits(cgObj1->getMassUnit(), cgObj1->getLengthUnit(),
		    cgObj1->getTimeUnit(), cgObj1->getTemperatureUnit(),
		    cgObj1->getAngleUnit());
  cgWrtObj->setBaseItrData(cgObj1->getBaseItrName(), 
                           cgObj1->getNTStep(), 
                           cgObj1->getTimeStep());
  cgWrtObj->setZoneItrData(cgObj1->getZoneItrName(), 
                           cgObj1->getGridCrdPntr(), 
                           cgObj1->getSolutionPntr());
  cgWrtObj->setZone(cgObj1->getZoneName(), cgObj1->getZoneType());
  cgWrtObj->setNVrtx(cgObj1->getNVertex());
  cgWrtObj->setNCell(cgObj1->getNElement());
  // define coordinates
  cgWrtObj->setGridXYZ(cgObj1->getVrtXCrd(), cgObj1->getVrtYCrd(), cgObj1->getVrtZCrd());
  // define connctivity
  cgWrtObj->setSection(cgObj1->getSectionName(), 
		       (ElementType_t) cgObj1->getElementType(), 
		       cgObj1->getElementConnectivity(-1));
  // define vertex and cell data 
  cgWrtObj->setSolutionNode("NodeData", Vertex);
  cgWrtObj->setSolutionNode("ElemData", CellCenter);

  // write skelleton of the file
  cgWrtObj->writeGridToFile();

  // writing data to the CGNS file
  std::vector<std::string> slnList;
  getSolutionDataNames(slnList);
  for (auto is=slnList.begin(); is!=slnList.end(); is++)
  {
    std::vector<double> trgSlnVec;
    solution_type_t st = trgCgObjs[0]->getSolutionData(*is, trgSlnVec);
    if (st== NODAL)      
      cgWrtObj->writeSolutionField(*is, "NodeData", RealDouble, &trgSlnVec[0]);
    else if (st == ELEMENTAL)
      cgWrtObj->writeSolutionField(*is, "ElemData", RealDouble, &trgSlnVec[0]);
  }
}

void gridTransfer::exportMeshToMAdLib(std::string gridName)
{
  // if gridName = source exports source mesh data to MAdLib, otherwise target
  // exporting mesh to the MAdLib
  if (!strcmp(gridName.c_str(), "src"))
  {
    MAd::GM_create(&srcModel,"");
    srcMesh = MAd::M_new(srcModel);
    exportToMAdMesh(srcMesh);
    classifyMAdMeshOpt(srcMesh);
    MAd::M_info(srcMesh, std::cout);
  } 
  else if (!strcmp(gridName.c_str(), "trg")) 
  {
    if (trgCgObjs.empty())
      return;
    MAd::GM_create(&trgModel,"");
    trgMesh = MAd::M_new(trgModel);
    trgCgObjs[0]->exportToMAdMesh(trgMesh);
    trgCgObjs[0]->classifyMAdMeshOpt(trgMesh);
    MAd::M_info(trgMesh, std::cout);
  } else {
    std::cerr << "Fatal Error: Only src or trg are accpeted.\n";
    throw;
  }
}

void gridTransfer::exportMeshBndToMAdLib(std::string gridName)
{
  // if gridName = source exports source mesh data to MAdLib, otherwise target
  // exporting mesh to the MAdLib
  if (!strcmp(gridName.c_str(), "src"))
  {
    MAd::GM_create(&srcModel,"");
    srcMesh = MAd::M_new(srcModel);
    exportToMAdMesh(srcMesh);
    classifyMAdMeshBnd(srcMesh);
    MAd::M_info(srcMesh, std::cout);
  } 
  else if (!strcmp(gridName.c_str(), "trg")) 
  {
    if (trgCgObjs.empty())
      return;
    MAd::GM_create(&trgModel,"");
    trgMesh = MAd::M_new(trgModel);
    trgCgObjs[0]->exportToMAdMesh(trgMesh);
    trgCgObjs[0]->classifyMAdMeshBnd(trgMesh);
    MAd::M_info(trgMesh, std::cout);
  } else {
    std::cerr << "Fatal Error: Only src or trg are accpeted.\n";
    throw;
  }
}

void gridTransfer::convertToMsh(std::string gridName)
{
  // if gridName = source exports source mesh data to MAdLib, otherwise target
  // exporting mesh to the MAdLib
  std::string fName;
  MAd::pMesh wrtMesh;

  if (!strcmp(gridName.c_str(), "src"))
  {
    if (!srcMesh) {std::cerr<<"Source needs to be exported to MAdLib format first.\n"; throw;}
    wrtMesh = srcMesh;
    fName = "source.msh";
  } 
  else if (!strcmp(gridName.c_str(), "trg")) 
  {
    if (!trgMesh) {std::cerr<<"Target needs to be exported to MAdLib format first.\n"; throw;}
    wrtMesh = trgMesh;
    fName = "target.msh";
  } 
  else 
  {
    std::cerr << "Fatal Error: Only src or trg are accpeted.\n";
    throw;
  }
  std::cout <<"Writing to gmsh format -> " << fName << std::endl;
  MAd::M_writeMsh(wrtMesh, fName.c_str(), 2, NULL);
  std::cout << __FILE__<<__LINE__<<std::endl;
}

void gridTransfer::convertToVTK(std::string gridName, bool withSolution)
{
  // if gridName = source exports source mesh data to MAdLib, otherwise target
  // exporting mesh to the MAdLib
  // first converting to msh format
  convertToMsh(gridName);
  // then convert grid to vtk
  std::string fName;
  MAd::pMesh wrtMesh;
  GModel* wrtGModel;
  wrtGModel = new GModel("default"); 

  if (!strcmp(gridName.c_str(), "src"))
  {
    wrtMesh = srcMesh;
    fName = "source.vtk";
    wrtGModel->readMSH("source.msh");
  } 
  else if (!strcmp(gridName.c_str(), "trg")) 
  {
    wrtMesh = trgMesh;
    fName = "target.vtk";
    wrtGModel->readMSH("target.msh");
  } 
  else 
  {
    std::cerr << "Fatal Error: Only src or trg are accpeted.\n";
    throw;
  }
  std::cout <<"Writing to vtk format -> " << fName << std::endl;
  wrtGModel->writeVTK(fName.c_str(), false, true);
  
  // writing solution data to the file
  if (withSolution)
  {
    if (!strcmp(gridName.c_str(), "src"))
    {
      vtkAnalyzer* wrtVTK;
      wrtVTK = new vtkAnalyzer((char *)fName.c_str());
      wrtVTK->read();
      // figure out what is existing on the stitched grid
      int outNData, outNDim;
      std::vector<std::string> slnNameList;
      std::vector<std::string> appSlnNameList;
      getSolutionDataNames(slnNameList);  
      getAppendedSolutionDataName(appSlnNameList);
      slnNameList.insert(slnNameList.end(),
			 appSlnNameList.begin(), appSlnNameList.end());
      // write all data into vtk file
      for (auto is=slnNameList.begin(); is<slnNameList.end(); is++)
      {
	std::vector<double> physData;
	getSolutionDataStitched(*is, physData, outNData, outNDim);
	solution_type_t dt = getSolutionDataObj(*is)->getDataType();
	if (dt == NODAL)      
	  wrtVTK->setPointDataArray((*is).c_str(), 1, physData);
	else if (dt == ELEMENTAL)
	  wrtVTK->setCellDataArray((*is).c_str(), 1, physData);
      }
      wrtVTK->report();
      wrtVTK->write((char *)fName.c_str());
    } 
    else if (!strcmp(gridName.c_str(), "trg"))
    {
      std::vector<std::string> appDataLst; 
      trgCgObjs[0]->getAppendedSolutionDataName(appDataLst);
      vtkAnalyzer* wrtVTK;
      wrtVTK = new vtkAnalyzer((char *)fName.c_str());
      wrtVTK->read();
      // write transfered data to vtk
      for (auto ia=appDataLst.begin(); ia!=appDataLst.end(); ia++)
      {
	std::vector<double> trgSlnVec;
	solution_type_t dt = trgCgObjs[0]->getSolutionData(*ia, trgSlnVec);
	if (dt == NODAL)
	  wrtVTK->setPointDataArray((*ia).c_str(), 1, trgSlnVec);
	else if (dt== ELEMENTAL)
	  wrtVTK->setCellDataArray((*ia).c_str(), 1, trgSlnVec);
	
      }
      wrtVTK->report();
      wrtVTK->write((char *)fName.c_str());
    }
  }
  
}

void gridTransfer::convertToSTL(std::string gridName, std::string prefix)
{
  // if gridName = src exports source mesh data to stl, otherwise trg
  // exports target mesh to stl
  // first converting to msh format
  //convertToMsh(gridName);

  // skin the gird
  std::string fName;
  std::cout << "Computing surface mesh.\n"; 
  std::vector<int> skinElmIds;
  MAd::pGModel tmpMdl = NULL;
  MAd::GM_create(&tmpMdl,"");
  MAd::pMesh skinMesh = M_new(tmpMdl);

  GModel* wrtGModel;
  wrtGModel = new GModel("default"); 

  if (!strcmp(gridName.c_str(), "src"))
  {
    srcMesh->skin_me(skinMesh, skinElmIds);
    skinMesh->classify_unclassified_entities();
    skinMesh->destroyStandAloneEntities();
    MAd::M_writeMsh(skinMesh, "srcSkinMesh.msh", 2, NULL);
    MAd::M_info(skinMesh);
    fName = "source.stl";
    wrtGModel->readMSH("srcSkinMesh.msh");
  } 
  else if (!strcmp(gridName.c_str(), "trg")) 
  {
    fName = "target.stl";
    wrtGModel->readMSH("target.msh");
  } 
  else 
  {
    std::cerr << "Fatal Error: Only src or trg are accpeted.\n";
    throw;
  }
  std::cout <<"Writing to stl format -> " << fName << std::endl;
  wrtGModel->writeSTL((prefix + fName).c_str(), false, true);
}

void gridTransfer::exportNodalDataToMAdLib()
{
  if (!srcMesh) {std::cerr<<"Source needs to be exported to MAdLib format first.\n"; throw;}
  MAd::NodalDataManagerSgl::instance().initialize(srcMesh); 
  // attach all nodal data
  int nData, nDim;
  std::vector<std::string> cgSlnNameList;
  std::vector<std::string> cgAppSlnNameList;
  getSolutionDataNames(cgSlnNameList);  
  getAppendedSolutionDataName(cgAppSlnNameList);
  cgSlnNameList.insert(cgSlnNameList.end(),
                     cgAppSlnNameList.begin(), cgAppSlnNameList.end());
  for (auto is=cgSlnNameList.begin(); is<cgSlnNameList.end(); is++)
  {
    std::vector<double> physData;
    getSolutionDataStitched(*is, physData, nData, nDim);
    solution_type_t dt = getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL) {
      std::cout << "MAdLib: Registering nodal data -> " << *is << std::endl;
      MAd::NodalDataManagerSgl::instance().registerData(*is, physData); 
    }
  }
  // report diagnostics
  MAd::NodalDataManagerSgl::instance().diagnostics(std::cout);
}

void gridTransfer::exportSrcToGModel()
{
  // TODO: Improve this method to avoid file I/O
  // convert to Msh format first
  convertToMsh("src");
  // read from file  
  srcGModel = new GModel("src");
  srcGModel->readMSH("source.msh");
  // some testing
  //std::cout << "Number of regions = " << srcGModel->getNumMeshElements() << std::endl; 
}

int gridTransfer::getElmIdx(std::vector<double>& xyz)
{
  // load to GModel if not yet
  if (!srcGModel)
    exportSrcToGModel();
  // find element indx containig the point
  SPoint3 pnt(xyz[0], xyz[1], xyz[2]);
  MElement* elm;
  elm = srcGModel->getMeshElementByCoord(pnt);
  if (!elm)
   return(-1);
  return(elm->getNum());
}

int gridTransfer::getBaryCrds(std::vector<double>& xyz, std::vector<double>& baryCrds, std::vector<int>& vrtIds)
{
  int elmIdx = getElmIdx(xyz);
  if (elmIdx == -1)
     return elmIdx;
  // get barycentric coords
  MAd::RIter rit = MAd::M_regionIter(srcMesh);
  MAd::pRegion pr;
  for (int iReg = 0; iReg<elmIdx; iReg++)
    pr = MAd::RIter_next(rit);
  double tmpBaryCrds[3];
  MAd::R_linearParams(pr, &(xyz[0]), tmpBaryCrds);
  baryCrds.push_back(1.0-tmpBaryCrds[0]-tmpBaryCrds[1]-tmpBaryCrds[2]);
  baryCrds.push_back(tmpBaryCrds[0]);
  baryCrds.push_back(tmpBaryCrds[1]);
  baryCrds.push_back(tmpBaryCrds[2]);
  // get vertex ids
  // vertex ids in MAdLib are start from 1
  // reduce one from it to make 0 indexed
  MAd::pPList rVerts = MAd::R_vertices(pr);
  void * tmp = NULL;
  while( MAd::pVertex pV = (MAd::pVertex)PList_next(rVerts,&tmp) )
    vrtIds.push_back(MAd::V_id(pV)-1); 
  /*
  // calculating point coordinates
  MAd::MDB_Point* pnt;
  double coords[3] = {0.0, 0.0, 0.0};
  double xyzv[4][3];
  R_coordP1(pr, xyzv);
  for (int i=0; i<4; i++)
    std::cout << xyzv[i][0] << " " << xyzv[i][1] << " " << xyzv[i][2] << std::endl;

  coords[0] = coords[0] + baryCrds[0]*xyzv[1][0];
  coords[1] = coords[1] + baryCrds[0]*xyzv[1][1];
  coords[2] = coords[2] + baryCrds[0]*xyzv[1][2];

  coords[0] = coords[0] + baryCrds[1]*xyzv[2][0];
  coords[1] = coords[1] + baryCrds[1]*xyzv[2][1];
  coords[2] = coords[2] + baryCrds[1]*xyzv[2][2];

  coords[0] = coords[0] + baryCrds[2]*xyzv[3][0];
  coords[1] = coords[1] + baryCrds[2]*xyzv[3][1];
  coords[2] = coords[2] + baryCrds[2]*xyzv[3][2];

  coords[0] = coords[0] + (1.0-baryCrds[0]-baryCrds[1]-baryCrds[2])*xyzv[0][0];
  coords[1] = coords[1] + (1.0-baryCrds[0]-baryCrds[1]-baryCrds[2])*xyzv[0][1];
  coords[2] = coords[2] + (1.0-baryCrds[0]-baryCrds[1]-baryCrds[2])*xyzv[0][2];
  std::cout << "Recons. Coords = "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<std::endl; 
  */
  return(elmIdx);
}

void gridTransfer::stitchMe(cgnsAnalyzer* cgObj, int zoneIdx, int verb)
{
  // load proper zone
  // asking object to load its zone information
  // this should refresh vertex coordinates and 
  // element connectivity tables
  cgObj->loadZone(zoneIdx, verb);

  // check if this is the first object being stitched
  if (nVertex==0)
  {
    // copying from the object
    isUnstructured = !cgObj->isStructured();
    cellDim = cgObj->getCellDim();
    physDim = cgObj->getPhysDim();
    massU = cgObj->getMassUnit();
    lengthU = cgObj->getLengthUnit();
    timeU = cgObj-> getTimeUnit();
    tempU = cgObj->getTemperatureUnit();
    angleU = cgObj->getAngleUnit();
    timeLabel = cgObj->getTimeStep();
    // using the very first zone as the initial
    // vertex and element informaiton
    nVertex = getZoneNVrtx(cgObj, 1);
    nElem = getZoneNCell(cgObj, 1);
    xCrd = getZoneCoords(cgObj, 1, 1);
    yCrd = getZoneCoords(cgObj, 1, 2);
    zCrd = getZoneCoords(cgObj, 1, 3);
    sectionType = (ElementType_t) getZoneRealSecType(cgObj, 1);
    elemConn = getZoneRealConn(cgObj, 1);
    stitchFldBc(cgObj, zoneIdx);
    return;
  }
  
  // (re)building the kdTree
  buildVertexKDTree();
  
  // clear old masks 
  vrtDataMask.clear();
  elmDataMask.clear();

  // adding new mesh non-repeating vertices
  std::vector<int>    newVrtIdx;
  std::vector<int>    rptVrtIdx;
  std::map<int,int>   rptVrtMap; // <newMeshIdx, currentMeshIdx>
  std::vector<double> newXCrd;
  std::vector<double> newYCrd;
  std::vector<double> newZCrd;
  int nNewVrt = 0;
  for (int iVrt=0; iVrt<cgObj->getNVertex(); iVrt++)
  {
    ANNpoint     qryVrtx;
    ANNidxArray  nnIdx;
    ANNdistArray dists;
    qryVrtx = annAllocPt(physDim);
    qryVrtx[0] = cgObj->getVrtXCrd(iVrt);
    qryVrtx[1] = cgObj->getVrtYCrd(iVrt);
    qryVrtx[2] = cgObj->getVrtZCrd(iVrt);
    nnIdx  = new ANNidx[1];
    dists  = new ANNdist[1];
    kdTree->annkSearch(qryVrtx, 1, nnIdx, dists);
    if (dists[0] > searchEps) {
      nNewVrt++;
      vrtDataMask.push_back(true);
      newVrtIdx.push_back(nVertex + nNewVrt);
      newXCrd.push_back(qryVrtx[0]);
      newYCrd.push_back(qryVrtx[1]);
      newZCrd.push_back(qryVrtx[2]);
    } else {
      vrtDataMask.push_back(false);
      newVrtIdx.push_back(nnIdx[0]+1);
      rptVrtIdx.push_back(iVrt);
      rptVrtMap[iVrt] = nnIdx[0]+1;
    }
  }
  std::cout << "Found " << nNewVrt << " new vertices.\n"; 
  std::cout << "Number of repeating index " << rptVrtIdx.size()
            << std::endl; 

  // currently implemented to add all new elements
  std::vector<int> newElemConn;
  int nNewElem=0;
  for (int iElem=0; iElem<cgObj->getNElement(); iElem++)
  {
    std::vector<int> rmtElemConn = cgObj->getElementConnectivity(iElem);
    
    // just adding all elements
    elmDataMask.push_back(true);
    nNewElem++;
    newElemConn.insert(newElemConn.end(), 
          rmtElemConn.begin(), rmtElemConn.end());
  }
  std::cout << "Found " << nNewElem << " new elements.\n";
   
  // switching conncetivity table to global
  for (int iIdx=0; iIdx<newElemConn.size(); iIdx++){
    newElemConn[iIdx] = newVrtIdx[newElemConn[iIdx]-1];
  }
  
  // stitching field values if requested
  stitchFldBc(cgObj, zoneIdx);

  // updating internal datastructure       
  nVertex += nNewVrt;
  nElem += nNewElem;
  xCrd.insert(xCrd.end(), newXCrd.begin(), newXCrd.end());
  yCrd.insert(yCrd.end(), newYCrd.begin(), newYCrd.end());
  zCrd.insert(zCrd.end(), newZCrd.begin(), newZCrd.end());
  elemConn.insert(elemConn.end(), newElemConn.begin(), newElemConn.end());
  zoneNames.push_back(cgObj->getFileName()+"->"+cgObj->getZoneName());
}

void gridTransfer::stitchFldBc(cgnsAnalyzer* cgObj, int zoneIdx, int verb)
{
  // load proper zone
  // asking object to load its zone information
  // this should refresh vertex coordinates and 
  // element connectivity tables
  cgObj->loadZone(zoneIdx);
  cgObj->clearAllSolutionData();
  // next two line just o triger  cgObj->populateSolutionDataNames();
  std::vector<std::string> listTmp;
  cgObj->getSolutionDataNames(listTmp);

  // treating special first time use case
  if (solutionMap.empty())
  {
    // populate solution data of current instance for 
    // stitching
    //cgObj->clearAllSolutionData();
    //cgObj->populateSolutionDataNames();
    solutionName = cgObj->getSolutionNodeNames();
    solutionGridLocation = cgObj->getSolutionGridLocations(); 
    solutionMap = cgObj->getSolutionMap();
    solutionNameLocMap = cgObj->getSolutionNameLocMap();
    // change current instance needed data and ask it
    // to load data into its container directly.
    indexFile = cgObj->getIndexFile();
    indexBase = cgObj->getIndexBase();
    indexZone = zoneIdx;
    loadSolutionDataContainer();

    // append Rocstar specific BCs as solution fields
    if (paneHasPatchNo(cgObj, zoneIdx))
    {
      appendSolutionData("patchNo", getPanePatchNo(cgObj, zoneIdx), 
			   ELEMENTAL, cgObj->getNElement(), 1); 
      appendSolutionData("bcflag", getPaneBcflag(cgObj, zoneIdx), 
			   ELEMENTAL, cgObj->getNElement(), 1); 
      appendSolutionData("cnstr_type", getPaneCnstrType(cgObj, zoneIdx), 
			   ELEMENTAL, cgObj->getNElement(), 1); 
    }
    return;
  }

  // Rocstar specific BCs remove the old exisiting field 
  if (paneHasPatchNo(cgObj, zoneIdx))
  {
    cgObj->delAppSlnData("patchNo");
    cgObj->appendSolutionData("patchNo", getPanePatchNo(cgObj, zoneIdx),
			       ELEMENTAL, cgObj->getNElement(), 1);
    cgObj->delAppSlnData("bcflag");
    cgObj->appendSolutionData("bcflag", getPaneBcflag(cgObj, zoneIdx),
			       ELEMENTAL, cgObj->getNElement(), 1);
    cgObj->delAppSlnData("cnstr_type");
    cgObj->appendSolutionData("cnstr_type", getPaneCnstrType(cgObj, zoneIdx),
			       ELEMENTAL, cgObj->getNElement(), 1);
  }

  // call current object stitch field
  stitchFields(cgObj);
}

void gridTransfer::stitchMe(gridTransfer* cgObj)
{ 
  std::cout << "Stitching rocstarCGNS to myself.\n"; 
  // check if this is the first object being stitched
  if (nVertex==0)
  {
    std::cerr << "This gridTransfer object needs to be initialized first.\n";
    exit(-1);
  }

  // (re)building the kdTree
  buildVertexKDTree();
  
  // clear old masks 
  vrtDataMask.clear();
  elmDataMask.clear();

  // adding new mesh non-repeating vertices
  std::vector<int>    newVrtIdx;
  std::vector<int>    rptVrtIdx;
  std::map<int,int>   rptVrtMap; // <newMeshIdx, currentMeshIdx>
  std::vector<double> newXCrd;
  std::vector<double> newYCrd;
  std::vector<double> newZCrd;
  int nNewVrt = 0;
  for (int iVrt=0; iVrt<cgObj->getNVertex(); iVrt++)
  {
    ANNpoint     qryVrtx;
    ANNidxArray  nnIdx;
    ANNdistArray dists;
    qryVrtx = annAllocPt(physDim);
    qryVrtx[0] = cgObj->getVrtXCrd(iVrt);
    qryVrtx[1] = cgObj->getVrtYCrd(iVrt);
    qryVrtx[2] = cgObj->getVrtZCrd(iVrt);
    nnIdx  = new ANNidx[1];
    dists  = new ANNdist[1];
    kdTree->annkSearch(qryVrtx, 1, nnIdx, dists);
    if (dists[0] > searchEps) {
      nNewVrt++;
      vrtDataMask.push_back(true);
      newVrtIdx.push_back(nVertex + nNewVrt);
      newXCrd.push_back(qryVrtx[0]);
      newYCrd.push_back(qryVrtx[1]);
      newZCrd.push_back(qryVrtx[2]);
    } else {
      vrtDataMask.push_back(false);
      newVrtIdx.push_back(nnIdx[0]+1);
      rptVrtIdx.push_back(iVrt);
      rptVrtMap[iVrt] = nnIdx[0]+1;
    }
  }
  std::cout << "Found " << nNewVrt << " new vertices.\n"; 
  std::cout << "Number of repeating index " << rptVrtIdx.size()
            << std::endl; 

  // currently implemented to add all new elements
  std::vector<int> newElemConn;
  int nNewElem=0;
  for (int iElem=0; iElem<cgObj->getNElement(); iElem++)
  {
    std::vector<int> rmtElemConn = cgObj->getElementConnectivity(iElem);
    
    // just adding all elements
    elmDataMask.push_back(true);
    nNewElem++;
    newElemConn.insert(newElemConn.end(), 
          rmtElemConn.begin(), rmtElemConn.end());
  }
  std::cout << "Found " << nNewElem << " new elements.\n";
   
  // switching conncetivity table to global
  for (int iIdx=0; iIdx<newElemConn.size(); iIdx++){
    newElemConn[iIdx] = newVrtIdx[newElemConn[iIdx]-1];
  }
  
  // stitching field values if requested
  //stitchFldBc(cgObj, zoneIdx);
  stitchFields((cgnsAnalyzer*) cgObj);

  // updating internal datastructure       
  nVertex += nNewVrt;
  nElem += nNewElem;
  xCrd.insert(xCrd.end(), newXCrd.begin(), newXCrd.end());
  yCrd.insert(yCrd.end(), newYCrd.begin(), newYCrd.end());
  zCrd.insert(zCrd.end(), newZCrd.begin(), newZCrd.end());
  elemConn.insert(elemConn.end(), newElemConn.begin(), newElemConn.end());
  //zoneNames.push_back(cgObj->getFileName()+"->"+cgObj->getZoneName());
  
}

////////////////////////////////////////////////////
// ZONE DATA ACCESS
////////////////////////////////////////////////////
int gridTransfer::getNZone(int indx)
{
  if (indx > nowCgFNames.size())
    return(-1);
  return(nowCgObjs[indx]->getNZone());
}

std::string gridTransfer::getZoneName(cgnsAnalyzer* cgObj, int zoneIdx)
{
  char zonename[33];
  int tmp[9];
  if (cg_zone_read(cgObj->getIndexFile(), 
               cgObj->getIndexBase(), 
               zoneIdx, zonename, tmp)) cg_error_exit();
  return(zonename);
}

std::string gridTransfer::getZoneName(int cgIdx, int zoneIdx)
{
  char zonename[33];
  int tmp[9];
  if (cg_zone_read(nowCgObjs[cgIdx]->getIndexFile(), 
               nowCgObjs[cgIdx]->getIndexBase(), 
               zoneIdx, zonename, tmp)) cg_error_exit();
  return(zonename);
}

ZoneType_t gridTransfer::getZoneType(int indx, int zidx)
{
  if (indx > nowCgFNames.size())
    return(ZoneTypeNull);
  nowCgObjs[indx]->loadZone(zidx);
  return(nowCgObjs[indx]->getZoneType());
}

std::string gridTransfer::getSectionName(int indx, int zidx)
{
  if (indx > nowCgFNames.size())
    return("INVALID");
  nowCgObjs[indx]->loadZone(zidx);
  return(nowCgObjs[indx]->getSectionName());
}

int gridTransfer::getElementType(int indx, int zidx)
{
  if (indx > nowCgFNames.size())
    return(-1);
  nowCgObjs[indx]->loadZone(zidx);
  return(nowCgObjs[indx]->getElementType());
}

int gridTransfer::getZoneNVrtx(cgnsAnalyzer* cgObj, int zoneIdx)
{
  char zonename[33];
  cgsize_t size[3];
  if (cg_zone_read(cgObj->getIndexFile(),
                   cgObj->getIndexBase(),
                   zoneIdx, zonename, size))
    cg_error_exit();
  return(size[0]);
}

int gridTransfer::getZoneNCell(cgnsAnalyzer* cgObj, int zoneIdx)
{
  char zonename[33];
  cgsize_t size[3];
  if (cg_zone_read(cgObj->getIndexFile(),
                   cgObj->getIndexBase(),
                   zoneIdx, zonename, size))
    cg_error_exit();
  return(size[1]);
}

std::vector<double> gridTransfer::getZoneCoords(cgnsAnalyzer* cgObj, int zoneIdx, int dim)
{
  std::vector<double> crds;
  if (cg_goto(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              "Zone_t", zoneIdx,
              "GridCoordinates_t", 1, "end"))
    cg_error_exit();
  char arrName[33];
  DataType_t dt;
  int dd;
  cgsize_t dimVec[3];
  if (cg_array_info(dim, arrName, &dt, &dd, dimVec)) cg_error_exit();
  crds.resize(dimVec[0],0);
  if (cg_array_read(dim, &crds[0])) cg_error_exit();
  // removing extra values that rocstar puts at the end
  for (int iEx=dimVec[0]; iEx>getZoneNVrtx(cgObj, zoneIdx); iEx--)
    crds.pop_back();
  return(crds);
}

std::vector<int> gridTransfer::getZoneRealConn(cgnsAnalyzer* cgObj, int zoneIdx)
{
  std::vector<int> conn;
  int secidx = 1;
  char secname[33];
  ElementType_t et;
  cgsize_t st, en;
  int nbndry, parflag;
  if (cg_section_read(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              zoneIdx, secidx,
              secname, &et, &st, &en, &nbndry, &parflag))
    cg_error_exit();
  // making sure we read real one otherwise next one is the real one
  if (std::string(secname).find("real") == std::string::npos)
    if (cg_section_read(cgObj->getIndexFile(),
	 cgObj->getIndexBase(),
	 zoneIdx, ++secidx,
	 secname, &et, &st, &en, &nbndry, &parflag))
      cg_error_exit();
  int nVrtxElm;
  switch(et)
  {
    case TETRA_4:
      nVrtxElm = 4;
      break;
    case HEXA_8:
      nVrtxElm = 8;
      break;
    case TRI_3:
      nVrtxElm = 3;
      break;
    case QUAD_4:
      nVrtxElm = 4;
      break;
    default:
      std::cerr << "Unknown element type " << st << std::endl;
      break;
  }
  conn.resize((en-st+1)*nVrtxElm,-1);
  if (cg_elements_read(cgObj->getIndexFile(), 
                       cgObj->getIndexBase(), 
                       zoneIdx, 
                       secidx, 
                       &conn[0], NULL))
    cg_error_exit();
  return(conn);  
}
 
int gridTransfer::getZoneRealSecType(cgnsAnalyzer* cgObj, int zoneIdx)
{
  std::vector<int> conn;
  int secidx = 1;
  char secname[33];
  ElementType_t et;
  cgsize_t st, en;
  int nbndry, parflag;
  if (cg_section_read(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              zoneIdx, secidx,
              secname, &et, &st, &en, &nbndry, &parflag))
    cg_error_exit();
  // making sure we read real one otherwise next one is the real one
  if (std::string(secname).find("real") == std::string::npos)
    if (cg_section_read(cgObj->getIndexFile(),
	 cgObj->getIndexBase(),
	 zoneIdx, ++secidx,
	 secname, &et, &st, &en, &nbndry, &parflag))
      cg_error_exit();
  return(et);
}

////////////////////////////////////////////////////
//  PANE DATA ACCESS (Example implementation of BCs support)
////////////////////////////////////////////////////
int gridTransfer::getPaneBcflag(cgnsAnalyzer* cgObj, int zoneIdx)
{
  if (cg_goto(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end")) 
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr=1; iArr<=nArr; iArr++)
  {
    char arrName[33];
    DataType_t dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "bcflag"))
      break;
  }
  if (iArr>nArr)
  {
    std::cerr << "Can not find bcflag." << std::endl;
    cg_error_exit();
  }
  int bcflag;
  if(cg_array_read(iArr, &bcflag)) cg_error_exit();
  return(bcflag);
}

bool gridTransfer::paneHasPatchNo(cgnsAnalyzer* cgObj, int zoneIdx)
{
  if (cg_goto(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end")) 
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr=1; iArr<=nArr; iArr++)
  {
    char arrName[33];
    DataType_t dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "patchNo"))
      return(true);
  }
  return(false);
}

int gridTransfer::getPanePatchNo(cgnsAnalyzer* cgObj, int zoneIdx)
{
  if (cg_goto(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end")) 
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr=1; iArr<=nArr; iArr++)
  {
    char arrName[33];
    DataType_t dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "patchNo"))
      break;
  }
  if (iArr>nArr)
  {
    std::cerr << "Can not find patchNo." << std::endl;
    cg_error_exit();
  }
  int patchNo;
  if(cg_array_read(iArr, &patchNo)) cg_error_exit();
  return(patchNo);
}

int gridTransfer::getPaneCnstrType(cgnsAnalyzer* cgObj, int zoneIdx)
{
  if (cg_goto(cgObj->getIndexFile(),
              cgObj->getIndexBase(),
              "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end")) 
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr=1; iArr<=nArr; iArr++)
  {
    char arrName[33];
    DataType_t dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "cnstr_type"))
      break;
  }
  if (iArr>nArr)
  {
    std::cerr << "Can not find patchNo." << std::endl;
    cg_error_exit();
  }
  int cnstrType;
  if(cg_array_read(iArr, &cnstrType)) cg_error_exit();
  return(cnstrType);
}

int gridTransfer::getNCgObj()
{
  return(nowCgObjs.size());
}

std::string gridTransfer::getBaseNameTrg()
{
  return(baseCgFNameTrg);
}

std::string gridTransfer::getBaseNameSrc()
{
  return(baseCgFNameSrc);
}

std::string gridTransfer::getBaseName(int indx)
{
  if (indx<nowCgObjs.size())
    return(nowCgObjs[indx]->getBaseName());
  return("INVALID");
}

std::string gridTransfer::getCgFName(int indx)
{
  if (indx < nowCgFNames.size())
    return(nowCgFNames[indx]);
  return("INVALID");
}

std::string gridTransfer::getBaseItrName(int indx)
{
  if (indx < nowCgFNames.size())
    return(nowCgObjs[indx]->getBaseItrName());
  return("INVALID");
}

int gridTransfer::getNTStep(int indx)
{
  if (indx < nowCgFNames.size())
    return(nowCgObjs[indx]->getNTStep());
  return(-1);
}

double gridTransfer::getTimeStep(int indx)
{
  if (indx < nowCgFNames.size())
    return(nowCgObjs[indx]->getTimeStep());
  return(-1);
}

std::string gridTransfer::getZoneItrName(int indx, int zidx)
{
  if (indx > nowCgFNames.size())
    return("INVALID");
  nowCgObjs[indx]->loadZone(zidx);
  return(nowCgObjs[indx]->getZoneItrName());
}

std::string gridTransfer::getGridCrdPntr(int indx, int zidx)
{
  if (indx > nowCgFNames.size())
    return("INVALID");
  nowCgObjs[indx]->loadZone(zidx);
  return(nowCgObjs[indx]->getGridCrdPntr());
}

std::string gridTransfer::getSolutionPntr(int indx, int zidx)
{
  if (indx > nowCgFNames.size())
    return("INVALID");
  nowCgObjs[indx]->loadZone(zidx);
  return(nowCgObjs[indx]->getSolutionPntr());
}
