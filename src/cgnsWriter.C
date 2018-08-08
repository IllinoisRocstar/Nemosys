/* Implementation of CGNSWriter class */
#include "cgnsWriter.H"
#include "GmshEntities.h"

void cgnsWriter::setUnits(MassUnits_t mu, LengthUnits_t lu, 
                          TimeUnits_t tu, TemperatureUnits_t tpu, AngleUnits_t au)
{
 massU = mu;
 lengthU = lu;
 timeU = tu;
 tempU = tpu;
 angleU = au;
}

void cgnsWriter::setBaseItrData(std::string bsitrname, int ntstp, double tval)
{
  baseItrName = bsitrname;
  nTStep = ntstp;
  timeLabel = tval;
}

void cgnsWriter::setIntData(std::string intname, int intval)
{
  intName = intname;
  intVal = intval;
}

void cgnsWriter::setZoneItrData(std::string zitrname, std::string grdptr, std::string slnptr)
{
  zoneItrName = zitrname;
  gridCrdPntr = grdptr;
  flowSlnPntr = slnptr;
}

void cgnsWriter::setZone(std::string zName, ZoneType_t zt)
{
  nZone++;
  zoneNames.push_back(zName);
  zoneName = zName;
  zoneTypes.push_back(zt);
  zoneType = zt; 
}

void cgnsWriter::setSection(std::string sName, ElementType_t st, vect<int>::v1d elmConn)
{
  nSection++;
  sectionNames.push_back(sName);
  sectionName = sName;
  sectionTypes.push_back(st);
  sectionType = st;
  elmConns.push_back(elmConn);
}

void cgnsWriter::setGlobalSection(std::string gsName, ElementType_t gst, vect<int>::v1d gelmConn)
{
  gnSection++;
  gsectionNames.push_back(gsName);
  gsectionName = gsName;
  gsectionTypes.push_back(gst);
  gsectionType = gst;
  gelmConns.push_back(gelmConn);
}

void cgnsWriter::setGlobalSection(std::string gsName, ElementType_t gst)
{
  gnSection++;
  gsectionNames.push_back(gsName);
  gsectionName = gsName;
  gsectionTypes.push_back(gst);
  gsectionType = gst;
  std::vector<int> tmp = {-1};
  gelmConns.push_back(tmp);
}

void cgnsWriter::resetSections()
{
  nSection = 0;
  sectionNames.clear();
  sectionName.clear();
  sectionTypes.clear();
  elmConns.clear();
  nCells.clear();
}

void cgnsWriter::resetGlobalSections()
{
  gnSection = 0;
  gsectionNames.clear();
  gsectionName.clear();
  gsectionTypes.clear();
  gelmConns.clear();
  gnCells.clear();
}

void cgnsWriter::setNVrtx(int nVt)
{
  nVrtx= nVt;
}

void cgnsWriter::setNCell(int nCl)
{
  nCell= nCl;
  nCells.push_back(nCl);
}

void cgnsWriter::setGlobalNCell(int gnCl)
{
  gnCell= gnCl;
  gnCells.push_back(gnCl);
}

void cgnsWriter::setGridXYZ(vect<double>::v1d x, vect<double>::v1d y, vect<double>::v1d z)
{
  xCrd.clear();
  xCrd.insert(xCrd.begin(), x.begin(), x.end());
  yCrd.clear();
  yCrd.insert(yCrd.begin(), y.begin(), y.end());
  zCrd.clear();
  zCrd.insert(zCrd.begin(), z.begin(), z.end());
}  

void cgnsWriter::setCoordRind(int rind)
{
  coordRind = rind;
}

void cgnsWriter::setVirtElmRind(int rind)
{
  virtElmRind = rind;
}

void cgnsWriter::setPconnGhostDescriptor(int ghostDescriptor)
{
  pConnGhostDescriptor = ghostDescriptor;
}

void cgnsWriter::setPconnVec(const vect<int>::v1d& _pConnVec)
{
  pConnVec = _pConnVec;
}

void cgnsWriter::setPatchNo(int _patchNo)
{
  patchNo = _patchNo;
  //std::cout << "patchNo = " << patchNo << std::endl;
}

void cgnsWriter::setVolCellFacesNumber(int _nVolCellFaces)
{
  nVolCellFaces = _nVolCellFaces;
}

void cgnsWriter::setBcflag(int _bcflag)
{
  bcflag = _bcflag;
  //std::cout << "bcflag = " << bcflag << std::endl;
}

void cgnsWriter::setCnstrtype(int _cnstr_type)
{
  cnstr_type = _cnstr_type;
  //std::cout << "cnstr_type = " << cnstr_type << std::endl;
}

void cgnsWriter::setSolutionNode(std::string ndeName, GridLocation_t slnLoc)
{
  if (slnLoc == Vertex)
    nVrtxSln++;
  else if (slnLoc == CellCenter)
    nCellSln++;
  else
    std::cerr << "Can not write to requested solution location.\n";
  solutionNameLocMap[ndeName] = slnLoc;
  slnNameNFld[ndeName] = 0;  
}

// write solution data node
void cgnsWriter::writeSolutionNode(std::string ndeName, GridLocation_t slnLoc)
{
  if (slnLoc == Vertex)
    nVrtxSln++;
  else if (slnLoc == CellCenter)
    nCellSln++;
  else
    std::cerr << "Can not write to requested solution location.\n";
  solutionNameLocMap[ndeName] = slnLoc;
  slnNameNFld[ndeName] = 0;  
  int slnIdx;
  if (cg_sol_write(indexFile, indexBase, indexZone, 
                   ndeName.c_str(), slnLoc, &slnIdx)) cg_error_exit();
  solutionIdx.push_back(slnIdx);
  solutionNameSolIdxMap[ndeName] = slnIdx;  
  cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", slnIdx, "end");
  cg_gridlocation_write(slnLoc);
 
  int rindArr[2]; 
  rindArr[0] = 0;
  if (coordRind && slnLoc == Vertex)
  {
    rindArr[1] = coordRind;
  }
  if (virtElmRind && slnLoc == CellCenter)
  {
    rindArr[1] = virtElmRind;
  }
  if (virtElmRind || coordRind)
  {
    if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,"FlowSolution_t",slnIdx,"end")) cg_error_exit();
    if (cg_rind_write(rindArr)) cg_error_exit();
  } 
}

/* Writing a solution field to the CGNS file.
   We assume the skeleton of the file is already written properly */
void cgnsWriter::writeSolutionField(std::string fname, std::string ndeName, DataType_t dt, void* data)
{
  std::cout << "writing out " << fname << std::endl;
  int slnIdx=-1;
  auto is = solutionNameLocMap.begin();
  while (is!=solutionNameLocMap.end())
  {
   slnIdx++;
   if (!strcmp((is->first).c_str(), ndeName.c_str()))
   {  
     break;
   }
   is++;
  }
  // sanity check
  if (is==solutionNameLocMap.end()){
    std::cerr << ndeName 
              << " is not an existing solution node.\n";
    return;
  }
  //std::cout << "1" << std::endl;
  // check vertex or cell based
  if (is->second == Vertex)
    nVrtxFld++;
  else
    nCellFld++;
  // write solution to file
  int fldIdx;
  //std::cout << "Sol Node Index : " << slnIdx << std::endl;
  //std::cout << "Node Name: " << ndeName << std::endl;
  //std::cout << "Sol name : " << fname << std::endl;
  //std::cout << "Index from name map : " << solutionNameSolIdxMap[ndeName] << std::endl;
  int currSlnIdx = solutionNameSolIdxMap[ndeName];
  if (cg_field_write(indexFile, indexBase, indexZone, currSlnIdx,//solutionIdx[slnIdx],
                     dt, fname.c_str(), data, &fldIdx)) cg_error_exit();
  //std::cout << "2" << std::endl;
  keyValueList fldIdxName;
  fldIdxName[fldIdx] = fname;
  std::pair<int, keyValueList> slnPair;
  slnPair.first = slnIdx;
  slnPair.second = fldIdxName;
  solutionMap[++nSlnFld] = slnPair; 
  // finding range of data
  double* tmpData = (double*) data;
  double min = tmpData[0];
  double max = tmpData[0];
  int nItr = 0;
  if (is->second == Vertex)
  {
    nItr = nVrtx;
  } else {
    nItr = nCell;
  }
  for (int it=0; it<nVrtx; it++) {
    min = std::min(tmpData[it], min);
    max = std::max(tmpData[it], max);
  }
  //std::cout << "3" << std::endl;
  // writing range descriptor
  std::ostringstream os;
  os << min << ", " << max;
  std::string range = os.str();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone,
              "FlowSolution_t", currSlnIdx, "DataArray_t", fldIdx, "end")) cg_error_exit();
  if (cg_descriptor_write("Range", range.c_str())) cg_error_exit();
  // write DimensionalExponents and units for cell data
  if (is->second == CellCenter)
  {
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", currSlnIdx,
                "DataArray_t", fldIdx, "end")) cg_error_exit();
    // dummy exponents 
    float exponents[5] = {0, 0, 0, 0, 0};
    if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
    if (cg_descriptor_write("Units", "dmy")) cg_error_exit();
  }
}

void cgnsWriter::writeGridToFile()
{
  // dummy variables
  std::ostringstream os;
  double min, max;
  std::string str;
  // writing base data
  char basename[33];
  strcpy(basename, baseName.c_str()); 
  if (cg_base_write(indexFile, basename, cellDim, physDim, &indexBase)) cg_error_exit();
  if (cg_goto(indexFile, indexBase,"end")) cg_error_exit();
  if (cg_units_write(massU, lengthU, timeU, tempU, angleU)) cg_error_exit();
  if (cg_biter_write(indexFile, indexBase, baseItrName.c_str(), nTStep)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, baseItrName.c_str(), 0, "end")) cg_error_exit();
  cgsize_t tmpDim = 1;
  if (cg_array_write("TimeValues", RealDouble, 1, (const cgsize_t*) 
                      &tmpDim, &timeLabel)) cg_error_exit();
  //writeZoneToFile();
}

void cgnsWriter::writeWinToFile()
{
  std::cout << "in writeWinToFile" << std::endl;
  intVal = 1;
  char basename[33];
  strcpy(basename, baseName.c_str()); 
  cgsize_t tmpDim[1] = {1};
  // create Integral data for window
  std::cout << "1" << std::endl;
  if (cg_goto(indexFile, indexBase, "end")) cg_error_exit();
  if (cg_integral_write(intName.c_str())) cg_error_exit();
  std::cout << "2" << std::endl;
  if (cg_goto(indexFile, indexBase, intName.c_str(), 0, "end")) cg_error_exit();
  double zoomFact_arr[1] = {intVal};
  if (cg_array_write("zoomFact", RealDouble, 1, tmpDim, zoomFact_arr)) cg_error_exit();
  std::cout << "3" << std::endl;
  if (cg_goto(indexFile, indexBase, intName.c_str(), 0,"zoomFact",0,"end")) cg_error_exit();
  std::ostringstream os;
  std::string str;
  os.str("");
  os.clear();
  os << std::to_string(intVal) << ", " << std::to_string(intVal);
  str = os.str();
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
  float dimUnits[5] = {0, 0, 0, 0, 0};
  if (cg_exponents_write(RealSingle, dimUnits)) cg_error_exit();
  if (cg_descriptor_write("Units", "none")) cg_error_exit();
  if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
}

//void cgnsWriter::writeVolGs(double* data)
//{
//  char basename[33];
//  strcpy(basename, baseName.c_str());
//  if (cg_base_read(indexFile, 1, basename, &cellDim, &physDim)) cg_error_exit();
//
//  double gsArray[nVolCellFaces] = {0};
//  cgsize_t tmpDim[1] = {nVolCellFaces};
//  if (cg_goto(indexFile, 1, "Zone_t", 1, "PaneData0", 0,"end")) cg_error_exit();
//  if (cg_array_write("gs", RealDouble, 1, tmpDim, data)) cg_error_exit();
//  if (cg_goto(indexFile, 1, "Zone_t", 1, "PaneData0",0,"gs",0,"end")) 
//    cg_error_exit();
//  if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
//  if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
//}

void cgnsWriter::writeZoneToFile(int typeFlag)
{
  // typeFlag
  // 0 = volume
  // 1 = rocflu surface
  // 2 = rocburn surface

  // dummy variables
  std::ostringstream os;
  double min, max;
  std::string str;
  // define zone name 
  char zonename[33];
  strcpy(zonename, zoneName.c_str());
  if (zoneType == Unstructured) {
    cgCoreSize[0]=nVrtx;
    cgCoreSize[1]=nCells[0];
    cgCoreSize[2]=0;
  } else {
    std::cerr << "Format is not supported.\n";
    cg_error_exit();
  }
  // create zone
  //std::cout << "ZONEe = " << zonename << std::endl;
  if (cg_zone_write(indexFile, indexBase, zonename, 
                    cgCoreSize, Unstructured, &indexZone)) cg_error_exit();

  // write zone iterative data
  if (cg_ziter_write(indexFile, indexBase, indexZone, zoneItrName.c_str())) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, 
                         zoneItrName.c_str(), 0, "end")) cg_error_exit();
  cgsize_t tmpDim2[2];
  tmpDim2[0] = 32;
  tmpDim2[1] = 1;
  if (cg_array_write("GridCoordinatesPointers", Character, 2,  tmpDim2, 
                     gridCrdPntr.c_str())) cg_error_exit();
  if (cg_array_write("FlowSolutionsPointers", Character, 2, tmpDim2, 
                     flowSlnPntr.c_str())) cg_error_exit();
 
  // create grid coordinates node
  int indexGrid = 1;
  if (cg_grid_write(indexFile, indexBase, indexZone, "GridCoordinates", &indexGrid)) cg_error_exit();
 
  // write rind
  if (coordRind)
  {
    int rindArr[2] = {0,coordRind};
    if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,"GridCoordinates",0,"end")) cg_error_exit();
    if (cg_rind_write(rindArr)) cg_error_exit();
  }

  // write coordinates + range, dimensional exponent, units
  float dimUnits[5] = {0, 1, 0, 0, 0};
  // x
  min = *std::min_element(xCrd.begin(), xCrd.end());
  max = *std::max_element(xCrd.begin(), xCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();
  if (cg_coord_write(indexFile,indexBase,indexZone,RealDouble,"CoordinateX",
       &xCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateX", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(RealSingle, dimUnits)) cg_error_exit();
  if (cg_descriptor_write("Units", "m")) cg_error_exit();
  // y
  min = *std::min_element(yCrd.begin(), yCrd.end());
  max = *std::max_element(yCrd.begin(), yCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();
  if (cg_coord_write(indexFile,indexBase,indexZone,RealDouble,"CoordinateY",
       &yCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateY", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(RealSingle, dimUnits)) cg_error_exit();
  if (cg_descriptor_write("Units", "m")) cg_error_exit();
  // z
  min = *std::min_element(zCrd.begin(), zCrd.end());
  max = *std::max_element(zCrd.begin(), zCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();
  if (cg_coord_write(indexFile,indexBase,indexZone,RealDouble,"CoordinateZ",
       &zCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateZ", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(RealSingle, dimUnits)) cg_error_exit();
  if (cg_descriptor_write("Units", "m")) cg_error_exit();

 
  // create link to the grid coordinates
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "end")) cg_error_exit();
  os.str("");
  os.clear();
  os <<"/"<<baseName<<"/"<<zoneName<<"/GridCoordinates";
  str = os.str();
  if (cg_link_write(gridCrdPntr.c_str(), "", str.c_str())) cg_error_exit();

  // write sections (real and virtual cells)
  int elm_start;
  for (int iSec=0; iSec<nSection; iSec++)
  {
    elm_start = (iSec == 0 ? 1 : elm_start+nCells[iSec-1]);
    int elm_end = elm_start + nCells[iSec] - 1;
    int nBdy = 0;
    if (sectionNames[iSec].compare("Empty:t3:virtual"))
    {
      int tmp_arr[3] = {0, 0, 0};
      if (cg_section_write(indexFile, indexBase, indexZone, (sectionNames[iSec]).c_str(),
                           sectionTypes[iSec], elm_start, elm_end, nBdy, 
                           &tmp_arr[0], &indexSection)) cg_error_exit();
    }
    else
    {
      if (cg_section_write(indexFile, indexBase, indexZone, (sectionNames[iSec]).c_str(),
                           sectionTypes[iSec], elm_start, elm_end, nBdy, 
                           &elmConns[iSec][0], &indexSection)) cg_error_exit();
    }
    
    std::cout << "on section " << iSec << " of " << nSection << "sections" << std::endl;
    // set virtual rind
    if (virtElmRind)
    {
      int rindArr[2] = {0,virtElmRind};

      if (!sectionNames[iSec].compare(":T4:virtual"))
      {
        if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,":T4:virtual",0,"end")) cg_error_exit();
        if (cg_rind_write(rindArr)) cg_error_exit();
      }
      else if (!sectionNames[iSec].compare(":T3:virtual"))
      {
        if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,":T3:virtual",0,"end")) cg_error_exit();
        if (cg_rind_write(rindArr)) cg_error_exit();
      }
    }
  }

  // create pane data node
  if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,0,"end")) cg_error_exit();
  if (cg_integral_write("PaneData0")) cg_error_exit();
  // write pconn (and ridges)
  if (!pConnVec.empty())
  {
    std::cout << "writing Pconn" << std::endl;
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0,"end")) cg_error_exit();
    cgsize_t tmpDim[1] = {pConnVec.size()};
    if (cg_array_write("pconn",Integer,1, tmpDim, &pConnVec[0])) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0, "pconn",0,"end")) 
      cg_error_exit();
    min = *std::min_element(pConnVec.begin(),pConnVec.end());
    max = *std::max_element(pConnVec.begin(),pConnVec.end());
    os.str("");
    os.clear();
    os << min << ", " << max;
    str = os.str();
    if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
    os.str("");
    os.clear();
    os << pConnGhostDescriptor;
    str = os.str();
    if (cg_descriptor_write("Ghost", str.c_str())) cg_error_exit();
    tmpDim[0] = 1;
    int ridge[1] = {0};
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0,"end")) cg_error_exit();
    if (cg_array_write("ridges#1of2", Integer,1,tmpDim,ridge)) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"ridges#1of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"end")) cg_error_exit();
    if (cg_array_write("ridges#2of2", Integer,1,tmpDim,ridge)) cg_error_exit();  
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"ridges#2of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit();

    // only for volume
    if (typeFlag == 0)
    {
      tmpDim[0] = nVolCellFaces;
      double gsArray[nVolCellFaces] = {0};
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"end")) cg_error_exit();
      if (cg_array_write("gs", RealDouble, 1, tmpDim, gsArray)) cg_error_exit();  
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"gs",0,"end")) 
        cg_error_exit();
      if (cg_descriptor_write("Range", "0, 0")) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
      // dummy exponents 
      float exponents[5] = {0, 0, 0, 0, 0};
      if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
      if (cg_descriptor_write("Units", "m/s")) cg_error_exit();
    }

    // only for surface
    if (typeFlag == 1)
    {
      std::cout << "writing bcflag" << std::endl;
      // bcflag
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0,"end")) cg_error_exit();
      int bcflag_arr[1] = {bcflag};
      if (cg_array_write("bcflag", Integer,1,tmpDim,bcflag_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"bcflag",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(bcflag) << ", " << std::to_string(bcflag);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
  
      std::cout << "writing patchNo" << std::endl;
      // patchNo
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0,"end")) cg_error_exit();
      int patchNo_arr[1] = {patchNo};
      if (cg_array_write("patchNo", Integer,1,tmpDim,patchNo_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"patchNo",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(patchNo) << ", " << std::to_string(patchNo);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
  
      std::cout << "writing cnstr_type" << std::endl;
      // cnstr_type
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0,"end")) cg_error_exit();
      int cnstr_type_arr[1] = {cnstr_type};
      if (cg_array_write("cnstr_type", Integer,1,tmpDim,cnstr_type_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,"cnstr_type",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(cnstr_type) << ", " << std::to_string(cnstr_type);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
  
      // write sections (global real and virtual surface cells)
      int elm_start;
      for (int igSec=0; igSec<gnSection; igSec++)
      {
        std::cout << "igSec = " << igSec << std::endl;
        elm_start = (igSec == 0 ? 1 : elm_start+gnCells[igSec-1]);
        std::cout << "elm_start = " << elm_start << std::endl;
        int elm_end;
        os.str("");
        os.clear();
        std::vector<int> elems;
        if (gnCells[igSec] == 0)
        {
          elm_end = elm_start;
          elems.push_back(0);
          os << "EMPTY";
        }
        else
        {
          elm_end = elm_start + gnCells[igSec] - 1;
          elems.push_back(1);
          os << std::to_string(min) << ", " << std::to_string(max);
        }
        str = os.str();
        int min, max;
        int elem_arr[elems.size()];
        int cnt = 0;
        for (auto it=elems.begin(); it != elems.end(); it++) {
          min = std::min(*it, min);
          max = std::max(*it, max);
          elem_arr[cnt] = *it;
          cnt++;
        }
        if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0", 0,"end")) cg_error_exit();
        if (cg_array_write(gsectionNames[igSec].c_str(), Integer, 1, tmpDim, elem_arr)) cg_error_exit();
        if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "PaneData0",0,gsectionNames[igSec].c_str(),0,"end")) 
          cg_error_exit();
        if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
        if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
      }
    }

  }
}


