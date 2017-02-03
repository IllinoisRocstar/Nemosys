/* Implementation of CGNSWriter class */
#include "cgnsWriter.H"

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

void cgnsWriter::setNVrtx(int nVt)
{
  nVrtx= nVt;
}

void cgnsWriter::setNCell(int nCl)
{
  nCell= nCl;
  nCells.push_back(nCl);
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

/* Writing a solution field to the CGNS file.
   We assume the skeleton of the file is already written properly */
void cgnsWriter::writeSolutionField(std::string fname, std::string ndeName, DataType_t dt, void* data)
{
  int slnIdx=-1;
  auto is = solutionNameLocMap.begin();
  while (is!=solutionNameLocMap.end())
  {
   slnIdx++;
   if (!strcmp((is->first).c_str(), ndeName.c_str()))
     break;
   is++;
  }
  // sanity check
  if (is==solutionNameLocMap.end()){
    std::cerr << ndeName 
              << " is not an existing solution node.\n";
    return;
  }
  // check vertex or cell based
  if (is->second == Vertex)
    nVrtxFld++;
  else
    nCellFld++;
  // write solution to file
  int fldIdx;
  if (cg_field_write(indexFile, indexBase, indexZone, solutionIdx[slnIdx],
                     dt, fname.c_str(), data, &fldIdx)) cg_error_exit();
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
  // writing range descriptor
  std::ostringstream os;
  os << min << ", " << max;
  std::string range = os.str();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone,
              "FlowSolution_t", solutionIdx[slnIdx], "DataArray_t", fldIdx, "end")) cg_error_exit();
  if (cg_descriptor_write("Range", range.c_str())) cg_error_exit();
  // write DimensionalExponents and units for cell data
  if (is->second == CellCenter)
  {
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", solutionIdx[slnIdx],
                "DataArray_t", fldIdx, "end")) cg_error_exit();
    // dummy exponents and units
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
  writeZoneToFile();
}

void cgnsWriter::writeZoneToFile()
{
  // dummy variables
  std::ostringstream os;
  double min, max;
  std::string str;
  // define zone name 
  char zonename[33];
  strcpy(zonename, zoneName.c_str());
  if (zoneType == Unstructured) {
    cgCoreSize[0]=nVrtx;
    cgCoreSize[1]=nCell;
    cgCoreSize[2]=0;
  } else {
    std::cerr << "Format is not supported.\n";
    cg_error_exit();
  }
  // create zone
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

  // write sections
  for (int iSec=0; iSec<nSection; iSec++)
  {
    int nBdy = 0;
    if (cg_section_write(indexFile, indexBase, indexZone, (sectionNames[iSec]).c_str(),
                         sectionTypes[iSec], 1, nCells[iSec], nBdy, 
                         &elmConns[iSec][0], &indexSection)) cg_error_exit();
  }
  // write solution data
  for (auto is=solutionNameLocMap.begin(); is!= solutionNameLocMap.end(); is++)
  {
    int slnIdx;
    if (cg_sol_write(indexFile, indexBase, indexZone, 
                     (is->first).c_str(), is->second, &slnIdx)) cg_error_exit();
    solutionIdx.push_back(slnIdx);
    cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", slnIdx, "end");
    cg_gridlocation_write(is->second);
  }

}


