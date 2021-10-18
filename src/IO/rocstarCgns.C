/* Special purpose class for Rocstar CGNS files */

#include "IO/rocstarCgns.H"

///////////////////////////////////////////////////
// INITIALIZATION
///////////////////////////////////////////////////

rocstarCgns::rocstarCgns(std::string fname)
    : cgnsAnalyzer(fname), myCgFName(fname), _burn(0) {
  std::size_t _loc = myCgFName.find_last_of("_");
  baseCgFName = myCgFName.substr(0, _loc + 1);
  padSize = myCgFName.size() - baseCgFName.size() - 5;
}

rocstarCgns::rocstarCgns(const std::vector<std::string> &fnames)
    : cgnsAnalyzer(fnames[0]), _burn(0), cgFNames(fnames) {
  std::size_t _loc = fnames[0].find_last_of("_");
  baseCgFName = fnames[0].substr(0, _loc + 1);
}

rocstarCgns::~rocstarCgns() {
  // cleaning up
  for (int ic = 0; ic < myCgObjs.size(); ic++) delete myCgObjs[ic];
}

/*
   Loads a series of cgns files based on current file name
   and rocstar's convention.
*/
void rocstarCgns::loadCgSeries(int nCg) {
  for (int iCg = 0; iCg < nCg; iCg++) {
    std::ostringstream suffix1;
    std::ostringstream suffix2;
    suffix1 << iCg;
    std::string tmp = suffix1.str();
    for (int iPad = 0; iPad < padSize - tmp.size(); iPad++) suffix2 << "0";
    suffix2 << tmp << ".cgns";
    std::string fName = baseCgFName + suffix2.str();
    cgnsAnalyzer *cgTmp = new cgnsAnalyzer(fName);
    cgTmp->loadGrid();
    myCgObjs.push_back(cgTmp);
    cgFNames.push_back(fName);
  }
}

void rocstarCgns::loadCgSeries() {
  for (int i = 0; i < cgFNames.size(); ++i) {
    cgnsAnalyzer *cgTmp = new cgnsAnalyzer(cgFNames[i]);
    cgTmp->loadGrid();
    myCgObjs.push_back(cgTmp);
  }
}

void rocstarCgns::closeCG() {
  for (auto it = myCgObjs.begin() + 1; it != myCgObjs.end(); it++)
    (*it)->closeCG();
}

int rocstarCgns::getNCgObj() { return (myCgObjs.size()); }

std::string rocstarCgns::getBaseName() { return (baseCgFName); }

std::string rocstarCgns::getBaseName(int indx) {
  if (indx < myCgObjs.size()) return (myCgObjs[indx]->getBaseName());
  return ("INVALID");
}

std::string rocstarCgns::getCgFName(int indx) {
  if (indx < cgFNames.size()) return (cgFNames[indx]);
  return ("INVALID");
}

std::string rocstarCgns::getBaseItrName(int indx) {
  if (indx < cgFNames.size()) return (myCgObjs[indx]->getBaseItrName());
  return ("INVALID");
}

int rocstarCgns::getNTStep(int indx) {
  if (indx < cgFNames.size()) return (myCgObjs[indx]->getNTStep());
  return (-1);
}

double rocstarCgns::getTimeStep(int indx) {
  if (indx < cgFNames.size()) return (myCgObjs[indx]->getTimeStep());
  return (-1);
}

std::string rocstarCgns::getZoneItrName(int indx, int zidx) {
  if (indx > cgFNames.size()) return ("INVALID");
  myCgObjs[indx]->loadZone(zidx);
  return (myCgObjs[indx]->getZoneItrName());
}

std::string rocstarCgns::getGridCrdPntr(int indx, int zidx) {
  if (indx > cgFNames.size()) return ("INVALID");
  myCgObjs[indx]->loadZone(zidx);
  return (myCgObjs[indx]->getGridCrdPntr());
}

std::string rocstarCgns::getSolutionPntr(int indx, int zidx) {
  if (indx > cgFNames.size()) return ("INVALID");
  myCgObjs[indx]->loadZone(zidx);
  return (myCgObjs[indx]->getSolutionPntr());
}
////////////////////////////////////////////////////
// DATA PROCESSING
////////////////////////////////////////////////////
void rocstarCgns::stitchGroup() {
  for (int iObj = 0; iObj < getNCgObj(); iObj++) {
    std::cout << "Reading " << myCgObjs[iObj]->getFileName() << std::endl;
    for (int iz = 1; iz <= myCgObjs[iObj]->getNZone(); iz++) {
      std::cout << "Zone " << getZoneName(myCgObjs[iObj], iz) << std::endl;
      stitchMe(myCgObjs[iObj], iz);
    }
  }
}

void rocstarCgns::stitchMe(cgnsAnalyzer *cgObj, int zoneIdx) {
  // load proper zone
  // asking object to load its zone information
  // this should refresh vertex coordinates and
  // element connectivity tables
  cgObj->loadZone(zoneIdx);

  // check if this is the first object being stitched
  if (nVertex == 0) {
    // copying from the object
    isUnstructured = !cgObj->isStructured();
    cellDim = cgObj->getCellDim();
    physDim = cgObj->getPhysDim();
    massU = cgObj->getMassUnit();
    lengthU = cgObj->getLengthUnit();
    timeU = cgObj->getTimeUnit();
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
    sectionType = (CGNS_ENUMT(ElementType_t))getZoneRealSecType(cgObj, 1);
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
  std::vector<int> newVrtIdx;
  std::vector<int> rptVrtIdx;
  std::map<int, int> rptVrtMap;  // <newMeshIdx, currentMeshIdx>
  std::vector<double> newXCrd;
  std::vector<double> newYCrd;
  std::vector<double> newZCrd;
  int nNewVrt = 0;
  for (int iVrt = 0; iVrt < cgObj->getNVertex(); iVrt++) {
    ANNpoint qryVrtx;
    ANNidxArray nnIdx;
    ANNdistArray dists;
    qryVrtx = annAllocPt(physDim);
    qryVrtx[0] = cgObj->getVrtXCrd(iVrt);
    qryVrtx[1] = cgObj->getVrtYCrd(iVrt);
    qryVrtx[2] = cgObj->getVrtZCrd(iVrt);
    nnIdx = new ANNidx[1];
    dists = new ANNdist[1];
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
      newVrtIdx.push_back(nnIdx[0] + 1);
      rptVrtIdx.push_back(iVrt);
      rptVrtMap[iVrt] = nnIdx[0] + 1;
    }
    delete[] nnIdx;
    delete[] dists;
    annDeallocPt(qryVrtx);
  }
  std::cout << "Found " << nNewVrt << " new vertices.\n";
  std::cout << "Number of repeating index " << rptVrtIdx.size() << std::endl;

  // currently implemented to add all new elements
  std::vector<int> newElemConn;
  int nNewElem = 0;
  for (int iElem = 0; iElem < cgObj->getNElement(); iElem++) {
    std::vector<cgsize_t> rmtElemConn = cgObj->getElementConnectivity(iElem);
    // just adding all elements
    elmDataMask.push_back(true);
    nNewElem++;
    newElemConn.insert(newElemConn.end(), rmtElemConn.begin(),
                       rmtElemConn.end());
  }
  std::cout << "Found " << nNewElem << " new elements.\n";

  // switching connectivity table to global
  for (int iIdx = 0; iIdx < newElemConn.size(); iIdx++) {
    newElemConn[iIdx] = newVrtIdx[newElemConn[iIdx] - 1];
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
  zoneNames.push_back(cgObj->getFileName() + "->" + cgObj->getZoneName());
}

void rocstarCgns::stitchFldBc(cgnsAnalyzer *cgObj, int zoneIdx) {
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
  if (solutionMap.empty()) {
    // populate solution data of current instance for
    // stitching
    // cgObj->clearAllSolutionData();
    // cgObj->populateSolutionDataNames();
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
    if (!_burn) {
      appendSolutionData("patchNo", getPanePatchNo(cgObj, zoneIdx),
                         solution_type_t::ELEMENTAL, cgObj->getNElement(), 1);
      appendSolutionData("bcflag", getPaneBcflag(cgObj, zoneIdx),
                         solution_type_t::ELEMENTAL, cgObj->getNElement(), 1);
      appendSolutionData("cnstr_type", getPaneCnstrType(cgObj, zoneIdx),
                         solution_type_t::ELEMENTAL, cgObj->getNElement(), 1);
    }
    return;
  }
  // Rocstar specific BCs remove the old exisiting field
  if (!_burn) {
    cgObj->delAppSlnData("patchNo");
    cgObj->appendSolutionData("patchNo", getPanePatchNo(cgObj, zoneIdx),
                              solution_type_t::ELEMENTAL, cgObj->getNElement(),
                              1);
    cgObj->delAppSlnData("bcflag");
    cgObj->appendSolutionData("bcflag", getPaneBcflag(cgObj, zoneIdx),
                              solution_type_t::ELEMENTAL, cgObj->getNElement(),
                              1);
    cgObj->delAppSlnData("cnstr_type");
    cgObj->appendSolutionData("cnstr_type", getPaneCnstrType(cgObj, zoneIdx),
                              solution_type_t::ELEMENTAL, cgObj->getNElement(),
                              1);
  }
  // call current object stitch field
  stitchFields(cgObj);
}

void rocstarCgns::stitchMe(rocstarCgns *cgObj) {
  std::cout << "Stitching rocstarCGNS to myself.\n";
  // check if this is the first object being stitched
  if (nVertex == 0) {
    std::cerr << "This rocstarCgns object needs to be initialized first.\n";
    exit(-1);
  }

  // (re)building the kdTree
  buildVertexKDTree();

  // clear old masks
  vrtDataMask.clear();
  elmDataMask.clear();

  // adding new mesh non-repeating vertices
  std::vector<int> newVrtIdx;
  std::vector<int> rptVrtIdx;
  std::map<int, int> rptVrtMap;  // <newMeshIdx, currentMeshIdx>
  std::vector<double> newXCrd;
  std::vector<double> newYCrd;
  std::vector<double> newZCrd;
  int nNewVrt = 0;
  for (int iVrt = 0; iVrt < cgObj->getNVertex(); iVrt++) {
    ANNpoint qryVrtx;
    ANNidxArray nnIdx;
    ANNdistArray dists;
    qryVrtx = annAllocPt(physDim);
    qryVrtx[0] = cgObj->getVrtXCrd(iVrt);
    qryVrtx[1] = cgObj->getVrtYCrd(iVrt);
    qryVrtx[2] = cgObj->getVrtZCrd(iVrt);
    nnIdx = new ANNidx[1];
    dists = new ANNdist[1];
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
      newVrtIdx.push_back(nnIdx[0] + 1);
      rptVrtIdx.push_back(iVrt);
      rptVrtMap[iVrt] = nnIdx[0] + 1;
    }
    delete[] nnIdx;
    delete[] dists;
    annDeallocPt(qryVrtx);
  }
  std::cout << "Found " << nNewVrt << " new vertices.\n";
  std::cout << "Number of repeating index " << rptVrtIdx.size() << std::endl;

  // currently implemented to add all new elements
  std::vector<int> newElemConn;
  int nNewElem = 0;
  for (int iElem = 0; iElem < cgObj->getNElement(); iElem++) {
    std::vector<cgsize_t> rmtElemConn = cgObj->getElementConnectivity(iElem);

    // just adding all elements
    elmDataMask.push_back(true);
    nNewElem++;
    newElemConn.insert(newElemConn.end(), rmtElemConn.begin(),
                       rmtElemConn.end());
  }
  std::cout << "Found " << nNewElem << " new elements.\n";

  // switching conncetivity table to global
  for (int iIdx = 0; iIdx < newElemConn.size(); iIdx++) {
    newElemConn[iIdx] = newVrtIdx[newElemConn[iIdx] - 1];
  }

  // stitching field values if requested
  // stitchFldBc(cgObj, zoneIdx);
  stitchFields((cgnsAnalyzer *)cgObj);

  // updating internal datastructure
  nVertex += nNewVrt;
  nElem += nNewElem;
  xCrd.insert(xCrd.end(), newXCrd.begin(), newXCrd.end());
  yCrd.insert(yCrd.end(), newYCrd.begin(), newYCrd.end());
  zCrd.insert(zCrd.end(), newZCrd.begin(), newZCrd.end());
  elemConn.insert(elemConn.end(), newElemConn.begin(), newElemConn.end());
  // zoneNames.push_back(cgObj->getFileName()+"->"+cgObj->getZoneName());
}

////////////////////////////////////////////////////
// ZONE DATA ACCESS
////////////////////////////////////////////////////
int rocstarCgns::getNZone(int indx) {
  if (indx > cgFNames.size()) return (-1);
  return (myCgObjs[indx]->getNZone());
}

std::string rocstarCgns::getZoneName(cgnsAnalyzer *cgObj, int zoneIdx) {
  char zonename[33];
  cgsize_t tmp[9];
  if (cg_zone_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                   zonename, tmp))
    cg_error_exit();
  return (zonename);
}

std::string rocstarCgns::getZoneName(int cgIdx, int zoneIdx) {
  char zonename[33];
  cgsize_t tmp[9];
  if (cg_zone_read(myCgObjs[cgIdx]->getIndexFile(),
                   myCgObjs[cgIdx]->getIndexBase(), zoneIdx, zonename, tmp))
    cg_error_exit();
  return (zonename);
}

CGNS_ENUMT(ZoneType_t) rocstarCgns::getZoneType(int indx, int zidx) {
  if (indx > cgFNames.size()) return (CGNS_ENUMV(ZoneTypeNull));
  myCgObjs[indx]->loadZone(zidx);
  return (myCgObjs[indx]->getZoneType());
}

std::string rocstarCgns::getSectionName(int indx, int zidx) {
  if (indx > cgFNames.size()) return ("INVALID");
  myCgObjs[indx]->loadZone(zidx);
  return (myCgObjs[indx]->getSectionName());
}

int rocstarCgns::getElementType(int indx, int zidx) {
  if (indx > cgFNames.size()) return (-1);
  myCgObjs[indx]->loadZone(zidx);
  return (myCgObjs[indx]->getElementType());
}

int rocstarCgns::getZoneNVrtx(cgnsAnalyzer *cgObj, int zoneIdx) {
  char zonename[33];
  cgsize_t size[3];
  if (cg_zone_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                   zonename, size))
    cg_error_exit();
  return (size[0]);
}

int rocstarCgns::getZoneNCell(cgnsAnalyzer *cgObj, int zoneIdx) {
  char zonename[33];
  cgsize_t size[3];
  if (cg_zone_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                   zonename, size))
    cg_error_exit();
  return (size[1]);
}

std::vector<double> rocstarCgns::getZoneCoords(cgnsAnalyzer *cgObj, int zoneIdx,
                                               int dim) {
  std::vector<double> crds;
  if (cg_goto(cgObj->getIndexFile(), cgObj->getIndexBase(), "Zone_t", zoneIdx,
              "GridCoordinates_t", 1, "end"))
    cg_error_exit();
  char arrName[33];
  CGNS_ENUMT(DataType_t) dt;
  int dd;
  cgsize_t dimVec[3];
  if (cg_array_info(dim, arrName, &dt, &dd, dimVec)) cg_error_exit();
  crds.resize(dimVec[0], 0);
  if (cg_array_read(dim, &crds[0])) cg_error_exit();
  // removing extra values that rocstar puts at the end
  for (int iEx = dimVec[0]; iEx > getZoneNVrtx(cgObj, zoneIdx); iEx--)
    crds.pop_back();
  return (crds);
}

std::vector<cgsize_t> rocstarCgns::getZoneRealConn(cgnsAnalyzer *cgObj,
                                                   int zoneIdx) {
  std::vector<cgsize_t> conn;
  int secidx = 1;
  char secname[33];
  CGNS_ENUMT(ElementType_t) et;
  cgsize_t st, en;
  int nbndry, parflag;
  if (cg_section_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                      secidx, secname, &et, &st, &en, &nbndry, &parflag))
    cg_error_exit();
  // making sure we read real one otherwise next one is the real one
  if (std::string(secname).find("real") == std::string::npos)
    if (cg_section_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                        ++secidx, secname, &et, &st, &en, &nbndry, &parflag))
      cg_error_exit();
  int nVrtxElm;
  switch (et) {
    case CGNS_ENUMV(TETRA_4): nVrtxElm = 4; break;
    case CGNS_ENUMV(HEXA_8): nVrtxElm = 8; break;
    case CGNS_ENUMV(TRI_3): nVrtxElm = 3; break;
    case CGNS_ENUMV(QUAD_4): nVrtxElm = 4; break;
    default: std::cerr << "Unknown element type " << st << std::endl; break;
  }
  conn.resize((en - st + 1) * nVrtxElm, -1);
  if (cg_elements_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                       secidx, &conn[0], NULL))
    cg_error_exit();
  return (conn);
}

int rocstarCgns::getZoneRealSecType(cgnsAnalyzer *cgObj, int zoneIdx) {
  std::vector<int> conn;
  int secidx = 1;
  char secname[33];
  CGNS_ENUMT(ElementType_t) et;
  cgsize_t st, en;
  int nbndry, parflag;
  if (cg_section_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                      secidx, secname, &et, &st, &en, &nbndry, &parflag))
    cg_error_exit();
  // making sure we read real one otherwise next one is the real one
  if (std::string(secname).find("real") == std::string::npos)
    if (cg_section_read(cgObj->getIndexFile(), cgObj->getIndexBase(), zoneIdx,
                        ++secidx, secname, &et, &st, &en, &nbndry, &parflag))
      cg_error_exit();
  return (et);
}

////////////////////////////////////////////////////
//  PANE DATA ACCESS
////////////////////////////////////////////////////
int rocstarCgns::getPaneBcflag(cgnsAnalyzer *cgObj, int zoneIdx) {
  if (cg_goto(cgObj->getIndexFile(), cgObj->getIndexBase(), "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end"))
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr = 1; iArr <= nArr; iArr++) {
    char arrName[33];
    CGNS_ENUMT(DataType_t) dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "bcflag")) break;
  }
  if (iArr > nArr) {
    std::cerr << "Can not find bcflag." << std::endl;
    cg_error_exit();
  }
  int bcflag;
  if (cg_array_read(iArr, &bcflag)) cg_error_exit();
  return (bcflag);
}

int rocstarCgns::getPanePatchNo(cgnsAnalyzer *cgObj, int zoneIdx) {
  if (cg_goto(cgObj->getIndexFile(), cgObj->getIndexBase(), "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end"))
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr = 1; iArr <= nArr; iArr++) {
    char arrName[33];
    CGNS_ENUMT(DataType_t) dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "patchNo")) break;
  }
  if (iArr > nArr) {
    std::cerr << "Can not find patchNo." << std::endl;
    cg_error_exit();
  }
  int patchNo;
  if (cg_array_read(iArr, &patchNo)) cg_error_exit();
  return (patchNo);
}

int rocstarCgns::getPaneCnstrType(cgnsAnalyzer *cgObj, int zoneIdx) {
  if (cg_goto(cgObj->getIndexFile(), cgObj->getIndexBase(), "Zone_t", zoneIdx,
              "IntegralData_t", 1, "end"))
    cg_error_exit();
  int nArr;
  int iArr;
  if (cg_narrays(&nArr)) cg_error_exit();
  for (iArr = 1; iArr <= nArr; iArr++) {
    char arrName[33];
    CGNS_ENUMT(DataType_t) dt;
    int dd;
    cgsize_t dimVec[3];
    if (cg_array_info(iArr, arrName, &dt, &dd, dimVec)) cg_error_exit();
    if (!strcmp(arrName, "cnstr_type")) break;
  }
  if (iArr > nArr) {
    std::cerr << "Can not find cnstr_type." << std::endl;
    cg_error_exit();
  }
  int cnstrType;
  if (cg_array_read(iArr, &cnstrType)) cg_error_exit();
  return (cnstrType);
}
