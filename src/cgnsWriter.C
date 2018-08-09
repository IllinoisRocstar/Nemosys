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

void cgnsWriter::setFluidUnitsMap()
{
  // fluid
  fluidUnitsMap.insert(std::pair<std::string,std::string>("CoordinateX",  "m"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("CoordinateY",  "m"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("CoordinateZ",  "m"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("dispX",        "none"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("dispY",        "none"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("dispZ",        "none"));
  //fluidUnitsMap.insert(std::pair<std::string,std::string>("pconn","none"));      
  //fluidUnitsMap.insert(std::pair<std::string,std::string>("ridges#1of2","none"));
  //fluidUnitsMap.insert(std::pair<std::string,std::string>("ridges#1of2","none"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("gs",           "m/s"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("rhof",         "kg/(m^3)"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("rhovfX",       "kg/(m^2 s)"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("rhovfY",       "kg/(m^2 s)"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("rhovfZ",       "kg/(m^2 s)"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("rhoEf",        "(J/(m^3))"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("pf",           "N/(m^2)"));
  fluidUnitsMap.insert(std::pair<std::string,std::string>("Tf",           "K"));
}

void cgnsWriter::setFluidDimMap()
{
  // fluid
  float dim1[5] = {0, 1, 0, 0, 0};
  float dim2[5] = {0, 0, 0, 0, 0};
  fluidDimMap.insert(std::pair<std::string,float*>("CoordinateX",   dim1));
  fluidDimMap.insert(std::pair<std::string,float*>("CoordinateY",   dim1));
  fluidDimMap.insert(std::pair<std::string,float*>("CoordinateZ",   dim1));
  //fluidDimMap.insert(std::pair<std::string,float*>("dispX",         dim));
  //fluidDimMap.insert(std::pair<std::string,float*>("dispY",         dim));
  //fluidDimMap.insert(std::pair<std::string,float*>("dispZ",         dim));
  //fluidDimMap.insert(std::pair<std::string,float*>("pconn", {}));      
  //fluidDimMap.insert(std::pair<std::string,float*>("ridges#1of2", {}));
  //fluidDimMap.insert(std::pair<std::string,float*>("ridges#1of2", {}));
  fluidDimMap.insert(std::pair<std::string,float*>("gs",            dim2));
  fluidDimMap.insert(std::pair<std::string,float*>("rhof",          dim2));
  fluidDimMap.insert(std::pair<std::string,float*>("rhovfX",        dim2));
  fluidDimMap.insert(std::pair<std::string,float*>("rhovfY",        dim2));
  fluidDimMap.insert(std::pair<std::string,float*>("rhovfZ",        dim2));
  fluidDimMap.insert(std::pair<std::string,float*>("rhoEf",         dim2)); 
  fluidDimMap.insert(std::pair<std::string,float*>("pf",            dim2));
  fluidDimMap.insert(std::pair<std::string,float*>("Tf",            dim2));
}

void cgnsWriter::setFluidMagMap()
{
  // fluid
  fluidMagMap.insert(std::pair<std::string,int>("CoordinateX",  0));
  fluidMagMap.insert(std::pair<std::string,int>("CoordinateY",  0));
  fluidMagMap.insert(std::pair<std::string,int>("CoordinateZ",  0));
  fluidMagMap.insert(std::pair<std::string,int>("dispX",        1));
  fluidMagMap.insert(std::pair<std::string,int>("dispY",        0));
  fluidMagMap.insert(std::pair<std::string,int>("dispZ",        0));
  //fluidMagMap.insert(std::pair<std::strinfloat*ng>("pconn", {}));
  //fluidMagMap.insert(std::pair<std::strinfloat*ng>("ridges#1of2",
  //fluidMagMap.insert(std::pair<std::strinfloat*ng>("ridges#1of2",
  fluidMagMap.insert(std::pair<std::string,int>("gs",           0));
  fluidMagMap.insert(std::pair<std::string,int>("rhof",         0));
  fluidMagMap.insert(std::pair<std::string,int>("rhovfX",       1));
  fluidMagMap.insert(std::pair<std::string,int>("rhovfY",       0));
  fluidMagMap.insert(std::pair<std::string,int>("rhovfZ",       0));
  fluidMagMap.insert(std::pair<std::string,int>("rhoEf",        0));
  fluidMagMap.insert(std::pair<std::string,int>("pf",           0));
  fluidMagMap.insert(std::pair<std::string,int>("Tf",           0));        
}

void cgnsWriter::setiFluidUnitsMap()
{
  // ifluid_ni/ifluid_b
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("CoordinateX", "m"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("CoordinateY", "m"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("CoordinateZ", "m"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("du_alpX",     "m"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("du_alpY",     "m"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("du_alpZ",     "m"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("vmX",         "m/s"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("vmY",         "m/s"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("vmZ",         "m/s"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("mdot_alp",    "kg/(m^2s)"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("rhofvf_alpX", "kg/(m^2s)"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("rhofvf_alpY", "kg/(m^2s)"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("rhofvf_alpZ", "kg/(m^2s)"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("Tflm_alp",    "K"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("Tb_alp",      "K"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("nf_alpX",     "none"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("nf_alpY",     "none"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("nf_alpZ",     "none"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("pf",          "Pa"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("tfX",         "Pa"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("tfY",         "Pa"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("tfZ",         "Pa"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("qc",          "kgK/(m^2s)"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("qr",          "kgK/(m^2s)"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("rhof_alp",    "kg/m^3"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("zoomFact",    "none"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("bflag",       "none"));
  ifluidUnitsMap.insert(std::pair<std::string,std::string>("mdot_old",    "kg/(m^2 s)"));
}


void cgnsWriter::setiFluidDimMap()
{
  // ifluid_ni/ifluid_b
  float dim1[5] = {0, 1, 0, 0, 0};
  float dim2[5] = {0, 0, 0, 0, 0};
  ifluidDimMap.insert(std::pair<std::string,float*>("CoordinateX",  dim1));
  ifluidDimMap.insert(std::pair<std::string,float*>("CoordinateY",  dim1));
  ifluidDimMap.insert(std::pair<std::string,float*>("CoordinateZ",  dim1));
  ifluidDimMap.insert(std::pair<std::string,float*>("du_alpX",      dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("du_alpY",      dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("du_alpZ",      dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("vmX",          dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("vmY",          dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("vmZ",          dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("mdot_alp",     dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("rhofvf_alpX",  dim2)); 
  ifluidDimMap.insert(std::pair<std::string,float*>("rhofvf_alpY",  dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("rhofvf_alpZ",  dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("Tflm_alp",     dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("Tb_alp",       dim2));
  //ifluidDimMap.insert(std::pair<std::string,float*>("nf_alpX",      dim));         
  //ifluidDimMap.insert(std::pair<std::string,float*>("nf_alpY",      dim));
  //ifluidDimMap.insert(std::pair<std::string,float*>("nf_alpZ",      dim));
  ifluidDimMap.insert(std::pair<std::string,float*>("pf",           dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("tfX",          dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("tfY",          dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("tfZ",          dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("qc",           dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("qr",           dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("rhof_alp",     dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("zoomFact",     dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("bflag",        dim2));
  ifluidDimMap.insert(std::pair<std::string,float*>("mdot_old",     dim2));
}

void cgnsWriter::setiFluidMagMap()
{
  // ifluid_ni/ifluid_b
  ifluidMagMap.insert(std::pair<std::string,int>("CoordinateX", 0));
  ifluidMagMap.insert(std::pair<std::string,int>("CoordinateY", 0));
  ifluidMagMap.insert(std::pair<std::string,int>("CoordinateZ", 0));
  ifluidMagMap.insert(std::pair<std::string,int>("du_alpX",     1));
  ifluidMagMap.insert(std::pair<std::string,int>("du_alpY",     0));
  ifluidMagMap.insert(std::pair<std::string,int>("du_alpZ",     0));
  ifluidMagMap.insert(std::pair<std::string,int>("vmX",         1));
  ifluidMagMap.insert(std::pair<std::string,int>("vmY",         0));
  ifluidMagMap.insert(std::pair<std::string,int>("vmZ",         0));
  ifluidMagMap.insert(std::pair<std::string,int>("mdot_alp",    0));
  ifluidMagMap.insert(std::pair<std::string,int>("rhofvf_alpX", 1));
  ifluidMagMap.insert(std::pair<std::string,int>("rhofvf_alpY", 0));
  ifluidMagMap.insert(std::pair<std::string,int>("rhofvf_alpZ", 0));
  ifluidMagMap.insert(std::pair<std::string,int>("Tflm_alp",    0));
  ifluidMagMap.insert(std::pair<std::string,int>("Tb_alp",      0));
  ifluidMagMap.insert(std::pair<std::string,int>("nf_alpX",     1));
  ifluidMagMap.insert(std::pair<std::string,int>("nf_alpY",     0));
  ifluidMagMap.insert(std::pair<std::string,int>("nf_alpZ",     0));
  ifluidMagMap.insert(std::pair<std::string,int>("pf",          0));
  ifluidMagMap.insert(std::pair<std::string,int>("tfX",         1));
  ifluidMagMap.insert(std::pair<std::string,int>("tfY",         0));
  ifluidMagMap.insert(std::pair<std::string,int>("tfZ",         0));
  ifluidMagMap.insert(std::pair<std::string,int>("qc",          0));
  ifluidMagMap.insert(std::pair<std::string,int>("qr",          0));
  ifluidMagMap.insert(std::pair<std::string,int>("rhof_alp",    0));
  ifluidMagMap.insert(std::pair<std::string,int>("zoomFact",    0));
  ifluidMagMap.insert(std::pair<std::string,int>("bflag",       0));
  ifluidMagMap.insert(std::pair<std::string,int>("mdot_old",    0));
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

void cgnsWriter::setTimestamp(std::string _trimmed_base_t)
{
  trimmed_base_t = _trimmed_base_t;
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
void cgnsWriter::writeSolutionNode(std::string ndeName, GridLocation_t slnLoc, int emptyFlag, int virtFlag, int typeFlag)
{
  // emptyFlag
  // 0 = write
  // 1 = don't write; node is empty

  // virtFlag
  // 0 = only write real
  // 1 = write real and virtual

  // typeFlag
  // 0 = surface
  // 1 = volume

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


  float* exponents;
  std::string units;
  if (typeFlag == 0)
  {
    if (ifluidDimMap.find("ndeName.c_str()") != ifluidDimMap.end())
    {
      exponents = ifluidDimMap[ndeName.c_str()];
      if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
    }
    if (ifluidDimMap.find("ndeName.c_str()") != ifluidDimMap.end())
    {
      units = ifluidUnitsMap[ndeName.c_str()];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
    }
  }
  else if (typeFlag == 1)
  {
    if (fluidDimMap.find("ndeName.c_str()") != fluidDimMap.end())
    {
      exponents = fluidDimMap[ndeName.c_str()];
      if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
    }
    if (fluidDimMap.find("ndeName.c_str()") != fluidDimMap.end())
    {
      units = fluidUnitsMap[ndeName.c_str()];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
    }
  }
  if (!emptyFlag)
  {
    solutionIdx.push_back(slnIdx);
    solutionNameSolIdxMap[ndeName] = slnIdx;  
    cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", slnIdx, "end");
    cg_gridlocation_write(slnLoc);
  
    if (virtFlag)
    {
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
  //std::cout << "sizeof(arr)/sizeof(*arr) = " << sizeof(data)/sizeof(*data) << std::endl;
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
  //if (is->second == CellCenter)
  //{
  //  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", currSlnIdx,
  //              "DataArray_t", fldIdx, "end")) cg_error_exit();
  //  // dummy exponents 
  //  float exponents[5] = {0, 0, 0, 0, 0};
  //  if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
  //  if (cg_descriptor_write("Units", "dmy")) cg_error_exit();
  //}
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
  if (cg_exponents_write(RealSingle, ifluidDimMap["zoomFact"])) cg_error_exit();
  if (cg_descriptor_write("Units", ifluidUnitsMap["zoomFact"].c_str())) cg_error_exit();
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
  std::cout << "writing zone type " << typeFlag << std::endl;

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

  float* exponents;
  std::string units;
  if (typeFlag == 0)
  {
    units = fluidUnitsMap["CoordinateX"];
    exponents = fluidDimMap["CoordinateX"];
  }
  else if (typeFlag == 1)
  {
    units = ifluidUnitsMap["CoordinateX"];
    exponents = ifluidDimMap["CoordinateX"];
  }
  std::cout << "22" << std::endl;
  std::cout << "exponents = " << exponents[0] << exponents[1] << exponents[2] << exponents[3] << exponents[4] << std::endl;
  if (cg_coord_write(indexFile,indexBase,indexZone,RealDouble,"CoordinateX",
       &xCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateX", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  std::cout << "2" << std::endl;
  std::cout << "exponents" << sizeof(exponents)/sizeof(*exponents) << std::endl;
  if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
  std::cout << "3" << std::endl;
  if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
  std::cout << "4" << std::endl;

  // y
  min = *std::min_element(yCrd.begin(), yCrd.end());
  max = *std::max_element(yCrd.begin(), yCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();
  if (typeFlag == 0)
  {
    units = fluidUnitsMap["CoordinateY"];
    exponents = fluidDimMap["CoordinateY"];
  }
  else if (typeFlag == 1)
  {
    units = ifluidUnitsMap["CoordinateY"];
    exponents = ifluidDimMap["CoordinateY"];
  }
  if (cg_coord_write(indexFile,indexBase,indexZone,RealDouble,"CoordinateY",
       &yCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateY", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
  if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
  // z
  min = *std::min_element(zCrd.begin(), zCrd.end());
  max = *std::max_element(zCrd.begin(), zCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();
  if (typeFlag == 0)
  {
    units = fluidUnitsMap["CoordinateZ"];
    exponents = fluidDimMap["CoordinateZ"];
  }
  else if (typeFlag == 1)
  {
    units = ifluidUnitsMap["CoordinateZ"];
    exponents = ifluidDimMap["CoordinateZ"];
  }
  if (cg_coord_write(indexFile,indexBase,indexZone,RealDouble,"CoordinateZ",
       &zCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateZ", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(RealSingle, exponents)) cg_error_exit();
  if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();

  std::cout << "3" << std::endl;

 
  // create link to the grid coordinates
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "end")) cg_error_exit();
  os.str("");
  os.clear();
  os <<"/"<<baseName<<"/"<<zoneName<<"/GridCoordinates";
  str = os.str();
  if (cg_link_write(gridCrdPntr.c_str(), "", str.c_str())) cg_error_exit();

  std::cout << "4" << std::endl;

  // write sections (real and virtual cells)
  int elm_start;
  for (int iSec=0; iSec<nSection; iSec++)
  {
    std::cout << "isec = " << iSec << std::endl;
    std::cout << "sectionNames[iSec]" << sectionNames[iSec] << std::endl;
    elm_start = (iSec == 0 ? 1 : elm_start+nCells[iSec-1]);
    int elm_end = elm_start + nCells[iSec] - 1;
    int nBdy = 0;
    std::cout << "comparison = " << !sectionNames[iSec].compare("Empty:t3:virtual") << std::endl;
    if (!sectionNames[iSec].compare("Empty:t3:virtual"))
    {
      std::cout << "5" << std::endl;
      int tmp_arr[3] = {0, 0, 0};
      if (cg_section_write(indexFile, indexBase, indexZone, (sectionNames[iSec]).c_str(),
                           sectionTypes[iSec], elm_start, elm_end, nBdy, 
                           &tmp_arr[0], &indexSection)) cg_error_exit();
    }
    else
    {
      std::cout << "6" << std::endl;
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
  std::string paneName("PaneData");
  paneName += this->trimmed_base_t;
  if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,0,"end")) cg_error_exit();
  if (cg_integral_write(paneName.c_str())) cg_error_exit();
  // write pconn (and ridges)
  cgsize_t tmpDim[1];
  if (!pConnVec.empty())
  {
    std::cout << "writing Pconn" << std::endl;
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
    tmpDim[0] = {pConnVec.size()};
    if (cg_array_write("pconn",Integer,1, tmpDim, &pConnVec[0])) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0, "pconn",0,"end")) 
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
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
    if (cg_array_write("ridges#1of2", Integer,1,tmpDim,ridge)) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"ridges#1of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"end")) cg_error_exit();
    if (cg_array_write("ridges#2of2", Integer,1,tmpDim,ridge)) cg_error_exit();  
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"ridges#2of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit();

    // only for volume
    if (typeFlag == 0)
    {
      tmpDim[0] = nVolCellFaces;
      double gsArray[nVolCellFaces] = {0};
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"end")) cg_error_exit();
      if (cg_array_write("gs", RealDouble, 1, tmpDim, gsArray)) cg_error_exit();  
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"gs",0,"end")) 
        cg_error_exit();
      if (cg_descriptor_write("Range", "0, 0")) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
      // dummy exponents 
      if (cg_exponents_write(RealSingle, fluidDimMap["gs"])) cg_error_exit();
      if (cg_descriptor_write("Units", fluidUnitsMap["gs"].c_str())) cg_error_exit();
    }

    // only for surface
    if (typeFlag == 1)
    {
      std::cout << "writing bcflag" << std::endl;
      // bcflag
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
      int bcflag_arr[1] = {bcflag};
      if (cg_array_write("bcflag", Integer,1,tmpDim,bcflag_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"bcflag",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(bcflag) << ", " << std::to_string(bcflag);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
  
      std::cout << "writing patchNo" << std::endl;
      // patchNo
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
      int patchNo_arr[1] = {patchNo};
      if (cg_array_write("patchNo", Integer,1,tmpDim,patchNo_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"patchNo",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(patchNo) << ", " << std::to_string(patchNo);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
  
      std::cout << "writing cnstr_type" << std::endl;
      // cnstr_type
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
      int cnstr_type_arr[1] = {cnstr_type};
      if (cg_array_write("cnstr_type", Integer,1,tmpDim,cnstr_type_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"cnstr_type",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(cnstr_type) << ", " << std::to_string(cnstr_type);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
  
      std::cout << "typeflag = " << typeFlag << std::endl;

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
        if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
        if (cg_array_write(gsectionNames[igSec].c_str(), Integer, 1, tmpDim, elem_arr)) cg_error_exit();
        if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,gsectionNames[igSec].c_str(),0,"end")) 
          cg_error_exit();
        if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
        if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
      }
    }

  }
  else if (typeFlag == 2)
  {
    std::cout << "writing pane data for rocburn" << std::endl;
    // write blank Pane Data
    tmpDim[0] = 1;
    int pconn[1] = {0};
    int ridge[1] = {0};

    std::cout << "1" << std::endl;

    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
    std::cout << "1" << std::endl;
    if (cg_array_write("pconn", Integer,1,tmpDim,ridge)) cg_error_exit();
    std::cout << "1" << std::endl;
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0, "pconn",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 

    std::cout << "2" << std::endl;


    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName, 0,"end")) cg_error_exit();
    if (cg_array_write("ridges#1of2", Integer,1,tmpDim,ridge)) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"ridges#1of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 

    std::cout << "3" << std::endl;

    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"end")) cg_error_exit();
    if (cg_array_write("ridges#2of2", Integer,1,tmpDim,ridge)) cg_error_exit();  
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName,0,"ridges#2of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
  }
}


