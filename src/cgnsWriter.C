/* Implementation of CGNSWriter class */
#include "cgnsWriter.H"
#include "GmshEntities.h"

// delete file
void cgnsWriter::deleteFile()
{
  std::remove(myCgFileName.c_str());
}

void cgnsWriter::setUnits(CG_MassUnits_t mu, CG_LengthUnits_t lu, 
                          CG_TimeUnits_t tu, CG_TemperatureUnits_t tpu, CG_AngleUnits_t au)
{
 massU = mu;
 lengthU = lu;
 timeU = tu;
 tempU = tpu;
 angleU = au;
}

// It is currently unclear whether or not Rocstar actually uses the coordinates or dimensional quantities.
// Currently they are just stored in the pre-defined maps below. The units defined below are copied directly
// from what is seen in Rocstar output files. Ideally, we could figure out a way to "transfer" this information
// from the old CGNS files to the new CGNS files. For now, we use the maps.

// Map defining fluid units
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
  fluidUnitsMap.insert(std::pair<std::string,std::string>("af",           "m/s"));
}

// Map defining fluid dimensions
void cgnsWriter::setFluidDimMap()
{
  // fluid
  std::vector<float> dim1 = {0, 1, 0, 0, 0};
  std::vector<float> dim2 = {0, 0, 0, 0, 0};
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateX",   dim1));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateY",   dim1));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateZ",   dim1));
  //fluidDimMap.insert(std::pair<std::string,float*>("dispX",         dim));
  //fluidDimMap.insert(std::pair<std::string,float*>("dispY",         dim));
  //fluidDimMap.insert(std::pair<std::string,float*>("dispZ",         dim));
  //fluidDimMap.insert(std::pair<std::string,float*>("pconn", {}));      
  //fluidDimMap.insert(std::pair<std::string,float*>("ridges#1of2", {}));
  //fluidDimMap.insert(std::pair<std::string,float*>("ridges#1of2", {}));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("gs",            dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhof",          dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhovfX",        dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhovfY",        dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhovfZ",        dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhoEf",         dim2)); 
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("pf",            dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("Tf",            dim2));
  fluidDimMap.insert(std::pair<std::string,std::vector<float>>("af",            dim2));
}

// Map storing a boolean that determines whether or not the "Magnitude" quantity is
// output to CGNS files. Currently, this is unused in the Remesh driver, but Rocstar
// outputs the magnitudes of some x-coordinate quantities, as seen below. This is not
// required functionality for Rocstar to run, however.
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
  fluidMagMap.insert(std::pair<std::string,int>("af",           0));              
}

// Map defining ifluid units
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

// Map defining ifluid dimensions
void cgnsWriter::setiFluidDimMap()
{
  // ifluid_ni/ifluid_b
  std::vector<float> dim1 = {0, 1, 0, 0, 0};
  std::vector<float> dim2 = {0, 0, 0, 0, 0};
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateX",  dim1));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateY",  dim1));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateZ",  dim1));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("du_alpX",      dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("du_alpY",      dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("du_alpZ",      dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("vmX",          dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("vmY",          dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("vmZ",          dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("mdot_alp",     dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhofvf_alpX",  dim2)); 
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhofvf_alpY",  dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhofvf_alpZ",  dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("Tflm_alp",     dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("Tb_alp",       dim2));
  //ifluidDimMap.insert(std::pair<std::string,float*>("nf_alpX",      dim));         
  //ifluidDimMap.insert(std::pair<std::string,float*>("nf_alpY",      dim));
  //ifluidDimMap.insert(std::pair<std::string,float*>("nf_alpZ",      dim));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("pf",           dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("tfX",          dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("tfY",          dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("tfZ",          dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("qc",           dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("qr",           dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("rhof_alp",     dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("zoomFact",     dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("bflag",        dim2));
  ifluidDimMap.insert(std::pair<std::string,std::vector<float>>("mdot_old",     dim2));
}

// Map storing magnitude quantities for ifluid files.
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

// Map defining burn units
void cgnsWriter::setBurnUnitsMap()
{
  // burn, iburn_all
  burnUnitsMap.insert(std::pair<std::string,std::string>("CoordinateX", "m"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("CoordinateY", "m"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("CoordinateZ", "m"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("bflag",       "none"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("pf",          "Pa"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("centersX",    "m"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("centersY",    "m"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("centersZ",    "m"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("rb",          "m/s"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("Tflm",        "K"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("rb_old",      "m/s"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("pf_old",      "Pa"));
  burnUnitsMap.insert(std::pair<std::string,std::string>("Tflm_old",        "K"));
}

// Map defining burn dimensions
void cgnsWriter::setBurnDimMap()
{
  // burn, iburn_all
  std::vector<float> dim1 = {0, 1, 0, 0, 0};
  std::vector<float> dim2 = {0, 0, 0, 0, 0};
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateX",  dim1));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateY",  dim1));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("CoordinateZ",  dim1));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("bflag",       dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("pf",          dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("centersX",    dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("centersY",    dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("centersZ",    dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("rb",          dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("Tflm",        dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("rb_old",      dim2)); 
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("pf_old",      dim2));
  burnDimMap.insert(std::pair<std::string,std::vector<float>>("Tflm_old",    dim2));
}

// Map storing magnitude quantities for burn files.
void cgnsWriter::setBurnMagMap()
{
  // burn, iburn_all
  burnMagMap.insert(std::pair<std::string,int>("CoordinateX", 0));
  burnMagMap.insert(std::pair<std::string,int>("CoordinateY", 0));
  burnMagMap.insert(std::pair<std::string,int>("CoordinateZ", 0));
  burnMagMap.insert(std::pair<std::string,int>("bflag",       0));
  burnMagMap.insert(std::pair<std::string,int>("pf",          0));
  burnMagMap.insert(std::pair<std::string,int>("centersX",    1));
  burnMagMap.insert(std::pair<std::string,int>("centersY",    0));
  burnMagMap.insert(std::pair<std::string,int>("centersZ",    0));
  burnMagMap.insert(std::pair<std::string,int>("rb",          0));
  burnMagMap.insert(std::pair<std::string,int>("Tflm",        0));
  burnMagMap.insert(std::pair<std::string,int>("rb_old",      0));
  burnMagMap.insert(std::pair<std::string,int>("pf_old",      0));
  burnMagMap.insert(std::pair<std::string,int>("Tflm_old",    0));
}


// Set type flag for CGNS writer
// 0 = volume (fluid)
// 1 = surface (ifluid)
// 2 = burn (iburn_all, burn)
void cgnsWriter::setTypeFlag(int _typeFlag)
{
  typeFlag = _typeFlag;
}

// Set base time iteration stamp string
void cgnsWriter::setBaseItrData(std::string bsitrname, int ntstp, double tval)
{
  baseItrName = bsitrname;
  nTStep = ntstp;
  timeLabel = tval;
}

// Set name and value for CGNS Integral writer
void cgnsWriter::setIntData(std::string intname, int intval)
{
  intName = intname;
  intVal = intval;
}

// Set CGNS Zone Iterative information
void cgnsWriter::setZoneItrData(std::string zitrname, std::string grdptr, std::string slnptr)
{
  zoneItrName = zitrname;
  gridCrdPntr = grdptr;
  flowSlnPntr = slnptr;
}

// Set zones to write to file
void cgnsWriter::setZone(std::string zName, CG_ZoneType_t zt)
{
  nZone++;
  zoneNames.push_back(zName);
  zoneName = zName;
  zoneTypes.push_back(zt);
  zoneType = zt; 
}

// Get number of sections
int cgnsWriter::getNSections()
{
  return sectionNames.size();
}

// Set local patch sections
void cgnsWriter::setSection(std::string sName, CG_ElementType_t st, vect<int>::v1d elmConn)
{
  nSection++;
  sectionNames.push_back(sName);
  sectionName = sName;
  sectionTypes.push_back(st);
  sectionType = st;
  elmConns.push_back(elmConn);
}

// Set global partition sections
void cgnsWriter::setGlobalSection(std::string gsName, CG_ElementType_t gst, vect<int>::v1d gelmConn)
{
  gnSection++;
  gsectionNames.push_back(gsName);
  gsectionName = gsName;
  gsectionTypes.push_back(gst);
  gsectionType = gst;
  gelmConns.push_back(gelmConn);
}

// Set global partition sections.
// This overload is for "blank" data sections. For this version of the driver which only supports
// tetrahedral meshes, these blank sections include hexahedral and quadrilateral element sections.
void cgnsWriter::setGlobalSection(std::string gsName, CG_ElementType_t gst)
{
  gnSection++;
  gsectionNames.push_back(gsName);
  gsectionName = gsName;
  gsectionTypes.push_back(gst);
  gsectionType = gst;
  std::vector<int> tmp = {0};
  gelmConns.push_back(tmp);
}

// Reset local patch sections
void cgnsWriter::resetSections()
{
  nSection = 0;
  sectionNames.clear();
  sectionName.clear();
  sectionTypes.clear();
  elmConns.clear();
  nCells.clear();
}

// Reset global partition sections
void cgnsWriter::resetGlobalSections()
{
  gnSection = 0;
  gsectionNames.clear();
  gsectionName.clear();
  gsectionTypes.clear();
  gelmConns.clear();
  gnCells.clear();
}

// Set timestamp to write to CGNS file
void cgnsWriter::setTimestamp(std::string _trimmed_base_t)
{
  trimmed_base_t = _trimmed_base_t;
}

// Set number of vertices to write to CGNS node data
void cgnsWriter::setNVrtx(int nVt)
{
  nVrtx= nVt;
}

// Set number of cells to write to CGNS element data
void cgnsWriter::setNCell(int nCl)
{
  nCell= nCl;
  nCells.push_back(nCl);
}

// Set global number of cells
void cgnsWriter::setGlobalNCell(int gnCl)
{
  gnCell= gnCl;
  gnCells.push_back(gnCl);
}

// Set grid coordinates, stored separately in x,y,z for CGNS
void cgnsWriter::setGridXYZ(vect<double>::v1d x, vect<double>::v1d y, vect<double>::v1d z)
{
  xCrd.clear();
  xCrd.insert(xCrd.begin(), x.begin(), x.end());
  yCrd.clear();
  yCrd.insert(yCrd.begin(), y.begin(), y.end());
  zCrd.clear();
  zCrd.insert(zCrd.begin(), z.begin(), z.end());
}  

// Set coordinate rind: virtual coordinates.
// The rule of thumb in CGNS files is: Total coordinates = number of vertexes + number of Rind
void cgnsWriter::setCoordRind(int rind)
{
  coordRind = rind;
}

// Set rind of virtual elements
void cgnsWriter::setVirtElmRind(int rind)
{
  virtElmRind = rind;
}

// Set number of ghost elements in Pane connectivity array
void cgnsWriter::setPconnGhostDescriptor(int ghostDescriptor)
{
  pConnGhostDescriptor = ghostDescriptor;
}

// Set entire Pane Connectivity vector
void cgnsWriter::setPconnVec(const vect<int>::v1d& _pConnVec)
{
  pConnVec = _pConnVec;
}

// Set limits, min and max, of Pane Connectivity vector, used in defining Pconn Range
void cgnsWriter::setPconnLimits(int _pconnProcMin, int _pconnProcMax)
{
  pConnMax = _pconnProcMax;
  pConnMin = _pconnProcMin;
}

// Set patch number
void cgnsWriter::setPatchNo(int _patchNo)
{
  patchNo = _patchNo;
}

// Set number of faces in the interior of the volume mesh, not including the boundary
// faces on the exterior of the volume mesh
void cgnsWriter::setVolCellFacesNumber(int _nVolCellFaces)
{
  nVolCellFaces = _nVolCellFaces;
}

// Set boundary condition flag
void cgnsWriter::setBcflag(int _bcflag)
{
  bcflag = _bcflag;
}

// Set cnstr type
void cgnsWriter::setCnstrtype(int _cnstr_type)
{
  cnstr_type = _cnstr_type;
}

// Sets names and locations of solutions to write
void cgnsWriter::setSolutionNode(std::string ndeName, CG_GridLocation_t slnLoc)
{
  if (slnLoc == CG_Vertex)
    nVrtxSln++;
  else if (slnLoc == CG_CellCenter)
    nCellSln++;
  else
    std::cerr << "Can not write to requested solution location.\n";
  solutionNameLocMap[ndeName] = slnLoc;
  slnNameNFld[ndeName] = 0;  
}

// write solution data node
void cgnsWriter::writeSolutionNode(std::string ndeName, CG_GridLocation_t slnLoc, int emptyFlag, int virtFlag)
{
  // emptyFlag
  // 0 = write
  // 1 = don't write; node is empty

  // virtFlag
  // 0 = only write real
  // 1 = write real and virtual

  if (slnLoc == CG_Vertex)
    nVrtxSln++;
  else if (slnLoc == CG_CellCenter)
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
  // surface (Rocflu fluid)
  if (typeFlag == 0)
  {
    if (fluidDimMap.find(ndeName) != fluidDimMap.end())
    {
      exponents = &fluidDimMap[ndeName.c_str()][0];
      if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
    }
    if (fluidDimMap.find(ndeName) != fluidDimMap.end())
    {
      units = fluidUnitsMap[ndeName.c_str()];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
    }
  }
  // volume (Rocflu ifluid)
  else if (typeFlag == 1)
  {
    if (ifluidDimMap.find(ndeName) != ifluidDimMap.end())
    {
      exponents = &ifluidDimMap[ndeName.c_str()][0];
      if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
    }
    if (ifluidDimMap.find(ndeName) != ifluidDimMap.end())
    {
      units = ifluidUnitsMap[ndeName.c_str()];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
    }
  }
  // burn (Rocburn iburn_all, burn)
  else if (typeFlag == 2)
  {
    if (burnDimMap.find(ndeName) != burnDimMap.end())
    {
      exponents = &burnDimMap[ndeName.c_str()][0];
      std::vector<float> test;
      test = burnDimMap[ndeName.c_str()];
      if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
    }
    if (burnDimMap.find(ndeName) != burnDimMap.end())
    {
      units = burnUnitsMap[ndeName.c_str()];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
    }
  }
  if (!emptyFlag)
  {
    solutionIdx.push_back(slnIdx);
    solutionNameSolIdxMap[ndeName] = slnIdx;  
    cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", slnIdx, "end");
    cg_gridlocation_write(slnLoc);
  
    // Write virtual
    if (virtFlag)
    {
      int rindArr[2]; 
      rindArr[0] = 0;
      if (coordRind && slnLoc == CG_Vertex)
      {
        rindArr[1] = coordRind;
      }
      if (virtElmRind && slnLoc == CG_CellCenter)
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
void cgnsWriter::writeSolutionField(std::string fname, std::string ndeName, CG_DataType_t dt, void* data)
{
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
  // check vertex or cell based
  if (is->second == CG_Vertex)
    nVrtxFld++;
  else
    nCellFld++;
  // write solution to file
  int fldIdx;
  int currSlnIdx = solutionNameSolIdxMap[ndeName];
  if (cg_field_write(indexFile, indexBase, indexZone, currSlnIdx,
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
  if (is->second == CG_Vertex)
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
              "FlowSolution_t", currSlnIdx, "DataArray_t", fldIdx, "end")) cg_error_exit();
  if (cg_descriptor_write("Range", range.c_str())) cg_error_exit();
  
  float* exponents;
  std::string units;
  // write DimensionalExponents and units for cell data
  if (is->second == CG_CellCenter)
  {
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", currSlnIdx,
                "DataArray_t", fldIdx, "end")) cg_error_exit();
    // surface (Rocflu fluid)
    if (typeFlag == 0)
    {
      if (fluidUnitsMap.find(fname) != fluidUnitsMap.end())
      {
        units = fluidUnitsMap[fname];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
      }
      if (fluidDimMap.find(fname) != fluidDimMap.end())
      {
        exponents = &fluidDimMap[fname][0];
        if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
      }
    }
    // volume (Rocflu ifluid)
    else if (typeFlag == 1)
    {
      if (ifluidUnitsMap.find(fname) != ifluidUnitsMap.end())
      {
        units = ifluidUnitsMap[fname];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
      }
      if (ifluidDimMap.find(fname) != ifluidDimMap.end())
      {
        exponents = &ifluidDimMap[fname][0];
        if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
      }
    }
    // burn (Rocburn iburn_all, burn)
    else if (typeFlag == 2)
    {
      if (burnUnitsMap.find(fname) != burnUnitsMap.end())
      {
        units = burnUnitsMap[fname];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
      }
      if (burnDimMap.find(fname) != burnDimMap.end())
      {
        exponents = &burnDimMap[fname][0];
        if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
      }
    }
  }
  // For vertex data only
  else if (is->second == CG_Vertex)
  {
    if (typeFlag == 1)
    {
      if (ifluidUnitsMap.find(fname) != ifluidUnitsMap.end())
      {
        units = ifluidUnitsMap[fname];
      if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
      }
      if (ifluidDimMap.find(fname) != ifluidDimMap.end())
      {
        exponents = &ifluidDimMap[fname][0];
        if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
      }
    }
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
  if (cg_array_write("TimeValues", CG_RealDouble, 1, (const cgsize_t*) 
                      &tmpDim, &timeLabel)) cg_error_exit();
}

void cgnsWriter::writeWinToFile()
{
  intVal = 1;
  char basename[33];
  strcpy(basename, baseName.c_str()); 
  cgsize_t tmpDim[1] = {1};
  // create Integral data for window
  if (cg_goto(indexFile, indexBase, "end")) cg_error_exit();
  if (cg_integral_write(intName.c_str())) cg_error_exit();
  if (cg_goto(indexFile, indexBase, intName.c_str(), 0, "end")) cg_error_exit();
  double zoomFact_arr[1] = {intVal};
  if (cg_array_write("zoomFact", CG_RealDouble, 1, tmpDim, zoomFact_arr)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, intName.c_str(), 0,"zoomFact",0,"end")) cg_error_exit();
  std::ostringstream os;
  std::string str;
  os.str("");
  os.clear();
  os << std::to_string(intVal) << ", " << std::to_string(intVal);
  str = os.str();
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
  float dimUnits[5] = {0, 0, 0, 0, 0};
  if (cg_exponents_write(CG_RealSingle, &ifluidDimMap["zoomFact"][0])) cg_error_exit();
  if (cg_descriptor_write("Units", ifluidUnitsMap["zoomFact"].c_str())) cg_error_exit();
  if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
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
  if (zoneType == CG_Unstructured) {
    cgCoreSize[0]=nVrtx;
    cgCoreSize[1]=nCells[0];
    cgCoreSize[2]=0;
  } else {
    std::cerr << "Format is not supported.\n";
    cg_error_exit();
  }

  // create zone
  if (cg_zone_write(indexFile, indexBase, zonename, 
                    cgCoreSize, CG_Unstructured, &indexZone)) cg_error_exit();

  // write zone iterative data
  if (cg_ziter_write(indexFile, indexBase, indexZone, zoneItrName.c_str())) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, 
                         zoneItrName.c_str(), 0, "end")) cg_error_exit();
  cgsize_t tmpDim2[2];
  tmpDim2[0] = 32;
  tmpDim2[1] = 1;
  if (cg_array_write("GridCoordinatesPointers", CG_Character, 2,  tmpDim2, 
                     gridCrdPntr.c_str())) cg_error_exit();
  if (cg_array_write("FlowSolutionsPointers", CG_Character, 2, tmpDim2, 
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

  float* exponents;
  std::string units;
  std::vector<float> test;

  // x-coordinates
  //---------------
  min = *std::min_element(xCrd.begin(), xCrd.end());
  max = *std::max_element(xCrd.begin(), xCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();

  // Get units and dimension information for volume, surface and burning fields
  if (typeFlag == 0)
  {
    units = fluidUnitsMap["CoordinateX"];
    exponents = &fluidDimMap["CoordinateX"][0];
    test = fluidDimMap["CoordinateX"];
  }
  else if (typeFlag == 1)
  {
    units = ifluidUnitsMap["CoordinateX"];
    exponents = &ifluidDimMap["CoordinateX"][0];
    test = ifluidDimMap["CoordinateX"];
  }
  else if (typeFlag == 2)
  {
    units = burnUnitsMap["CoordinateX"];
    exponents = &burnDimMap["CoordinateX"][0];
    test = burnDimMap["CoordinateX"];
  }
  if (cg_coord_write(indexFile,indexBase,indexZone,CG_RealDouble,"CoordinateX",
       &xCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateX", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
  if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();

  // y-coordinates
  //---------------
  min = *std::min_element(yCrd.begin(), yCrd.end());
  max = *std::max_element(yCrd.begin(), yCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();

  // Get units and dimension information for volume, surface and burning fields
  if (typeFlag == 0)
  {
    units = fluidUnitsMap["CoordinateY"];
    exponents = &fluidDimMap["CoordinateY"][0];
  }
  else if (typeFlag == 1)
  {
    units = ifluidUnitsMap["CoordinateY"];
    exponents = &ifluidDimMap["CoordinateY"][0];
  }
  else if (typeFlag == 2)
  {
    units = burnUnitsMap["CoordinateY"];
    exponents = &burnDimMap["CoordinateY"][0];
  }
  if (cg_coord_write(indexFile,indexBase,indexZone,CG_RealDouble,"CoordinateY",
       &yCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateY", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
  if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();

  // z-coordinates
  //---------------
  min = *std::min_element(zCrd.begin(), zCrd.end());
  max = *std::max_element(zCrd.begin(), zCrd.end());
  os.str("");
  os.clear();
  os << min << ", " << max;
  str = os.str();

  // Get units and dimension information for volume, surface and burning fields
  if (typeFlag == 0)
  {
    units = fluidUnitsMap["CoordinateZ"];
    exponents = &fluidDimMap["CoordinateZ"][0];
  }
  else if (typeFlag == 1)
  {
    units = ifluidUnitsMap["CoordinateZ"];
    exponents = &ifluidDimMap["CoordinateZ"][0];
  }
  else if (typeFlag == 2)
  {
    units = burnUnitsMap["CoordinateZ"];
    exponents = &burnDimMap["CoordinateZ"][0];
  }
  if (cg_coord_write(indexFile,indexBase,indexZone,CG_RealDouble,"CoordinateZ",
       &zCrd[0],&indexCoord)) cg_error_exit();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "GridCoordinates", 0,
              "CoordinateZ", 0, "end")) cg_error_exit(); 
  if (cg_descriptor_write("Range", str.c_str())) cg_error_exit();
  if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
  if (cg_descriptor_write("Units", units.c_str())) cg_error_exit();
 
  // create link to the grid coordinates
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "end")) cg_error_exit();
  os.str("");
  os.clear();
  os <<"/"<<baseName<<"/"<<zoneName<<"/GridCoordinates";
  str = os.str();
  if (cg_link_write(gridCrdPntr.c_str(), "", str.c_str())) cg_error_exit();

  // write Sections (real and virtual cells)
  int elm_start;
  int elm_end;
  int nBdy;
  for (int iSec=0; iSec<nSection; iSec++)
  {
    elm_start = (iSec == 0 ? 1 : elm_start+nCells[iSec-1]);
    int elm_end = elm_start + nCells[iSec] - 1;
    int nBdy = 0;
    if (!sectionNames[iSec].compare("Empty:t3:virtual"))
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
    // Write virtual elements
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
      else if (!sectionNames[iSec].compare(":t3:virtual"))
      {
        if (cg_goto(indexFile, indexBase, "Zone_t",indexZone,":t3:virtual",0,"end")) cg_error_exit();
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
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
    tmpDim[0] = {pConnVec.size()};
    if (cg_array_write("pconn",CG_Integer,1, tmpDim, &pConnVec[0])) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0, "pconn",0,"end")) 
      cg_error_exit();
    min = pConnMin;
    max = pConnMax;

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
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
    if (cg_array_write("ridges#1of2", CG_Integer,1,tmpDim,ridge)) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"ridges#1of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"end")) cg_error_exit();
    if (cg_array_write("ridges#2of2", CG_Integer,1,tmpDim,ridge)) cg_error_exit();  
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"ridges#2of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit();

    // only for volume
    if (typeFlag == 0)
    {
      tmpDim[0] = nVolCellFaces;
      double gsArray[nVolCellFaces] = {0};
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"end")) cg_error_exit();
      if (cg_array_write("gs", CG_RealDouble, 1, tmpDim, gsArray)) cg_error_exit();  
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"gs",0,"end")) 
        cg_error_exit();
      if (cg_descriptor_write("Range", "0, 0")) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
      // dummy exponents 
      if (cg_exponents_write(CG_RealSingle, &fluidDimMap["gs"][0])) cg_error_exit();
      if (cg_descriptor_write("Units", fluidUnitsMap["gs"].c_str())) cg_error_exit();
    }

    // only for surface
    if (typeFlag == 1)
    {
      // bcflag
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
      int bcflag_arr[1] = {bcflag};
      if (cg_array_write("bcflag", CG_Integer,1,tmpDim,bcflag_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"bcflag",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(bcflag) << ", " << std::to_string(bcflag);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
  
      // patchNo
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
      int patchNo_arr[1] = {patchNo};
      if (cg_array_write("patchNo", CG_Integer,1,tmpDim,patchNo_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"patchNo",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(patchNo) << ", " << std::to_string(patchNo);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 
  
      // cnstr_type
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
      int cnstr_type_arr[1] = {cnstr_type};
      if (cg_array_write("cnstr_type", CG_Integer,1,tmpDim,cnstr_type_arr)) cg_error_exit();
      if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"cnstr_type",0,"end")) 
        cg_error_exit();
      os.str("");
      os.clear();
      os << std::to_string(cnstr_type) << ", " << std::to_string(cnstr_type);
      str = os.str();
      if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 
      if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
    }

  }
  // only for burn
  else if (typeFlag == 2)
  {
    // write blank Pane Data
    tmpDim[0] = 1;
    int pconn[1] = {0};
    int ridge[1] = {0};

    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
    if (cg_array_write("pconn", CG_Integer,1,tmpDim,ridge)) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0, "pconn",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 

    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
    if (cg_array_write("ridges#1of2", CG_Integer,1,tmpDim,ridge)) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"ridges#1of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit(); 

    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"end")) cg_error_exit();
    if (cg_array_write("ridges#2of2", CG_Integer,1,tmpDim,ridge)) cg_error_exit();  
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(),0,"ridges#2of2",0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", "EMPTY")) cg_error_exit(); 
    if (cg_descriptor_write("Ghost", "0")) cg_error_exit();
  }

  // Write global sections. Global here means indexing according to partition, instead of local to
  // patch (as above).
  // There are far too many conditional statements below. This could be made neater, but the output
  // format of the global surface sections is a little bit weird. It would be worth exploring how
  // much of this "weirdness" Rocstar actually depends on.
  for (int igSec=0; igSec<gnSection; igSec++)
  {
    elm_start = 1;
    if (gnCells[igSec] == 0)
    {
      elm_end = elm_start;
    }
    else
    {
      elm_end = gnCells[igSec];
    }
    nBdy = 0;

    if (gsectionNames[igSec] == "t3g:virtual#1of3" || 
        gsectionNames[igSec] == "t3g:virtual#2of3" ||
        gsectionNames[igSec] == "t3g:virtual#3of3")
    {
      min = gelmConns[igSec][0];
      max = gelmConns[igSec][0];
      os.str("");
      os.clear();
      if (min == 0 && max == 0)
      {
        os << "EMPTY";
      }
      else
      {
        os << min << ", " << max;
      }
      str = os.str();
    }
    else if (gsectionNames[igSec] == "q4g:real#1of4" || 
             gsectionNames[igSec] == "q4g:real#2of4" ||
             gsectionNames[igSec] == "q4g:real#3of4" ||
             gsectionNames[igSec] == "q4g:real#4of4" ||
             gsectionNames[igSec] == "q4g:virtual#1of4" ||
             gsectionNames[igSec] == "q4g:virtual#2of4" ||
             gsectionNames[igSec] == "q4g:virtual#3of4" ||
             gsectionNames[igSec] == "q4g:virtual#4of4")
    {
      os.str("");
      os.clear();
      os << "EMPTY";
      str = os.str();
    }
    else
    {
      min = *std::min_element(gelmConns[igSec].begin(),gelmConns[igSec].end());
      max = *std::max_element(gelmConns[igSec].begin(),gelmConns[igSec].end());
      os.str("");
      os.clear();
      os << min << ", " << max;
      str = os.str();
    }

    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0,"end")) cg_error_exit();
    tmpDim[0] = {gelmConns[igSec].size()};
    if (cg_array_write((gsectionNames[igSec]).c_str(),CG_Integer,1, tmpDim, &gelmConns[igSec][0])) cg_error_exit();
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, paneName.c_str(), 0, (gsectionNames[igSec]).c_str(),0,"end")) 
      cg_error_exit();
    if (cg_descriptor_write("Range", str.c_str())) cg_error_exit(); 

    if (gsectionNames[igSec] == "t3g:virtual#1of3" || 
        gsectionNames[igSec] == "t3g:virtual#2of3" ||
        gsectionNames[igSec] == "t3g:virtual#3of3")
    {
      os.str("");
      os.clear();
      if (gnCells[igSec] == 1 && (min == 0 && max == 0))
      {
        os << 0;
      }
      else
      {
        os << elm_end;
      }
      str = os.str();
    }
    else
    {
      os.str("");
      os.clear();
      os << 0;
      str = os.str();
    }
    if (cg_descriptor_write("Ghost", str.c_str() )) cg_error_exit();
  }
}
