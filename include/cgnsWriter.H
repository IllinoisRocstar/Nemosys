/* a class for writing CGNS files */
#ifndef __CGNSWRITER_H
#define __CGNSWRITER_H

// others
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream> 

//MAd Lib
#include "ModelInterface.h"
#include "MeshDataBaseInterface.h"
//#include "GmshEntities.h"

// third party
#include "cgnslib.h"

// Nemosys
#include "cgnsAnalyzer.H"

// typedefs
template<class T>
struct vect{
  typedef std::vector<T> v1d;
  typedef std::vector<std::vector<T> > v2d;
};

class cgnsWriter{
public:
  cgnsWriter(std::string fname, std::string baseName, int cellDim, int physDim) :
  myCgFileName(fname), cellDim(cellDim), physDim(physDim),
  baseName(baseName), nBase(1), nZone(0), nSection(0),
  nVrtx(0), nCell(0), nSlnFld(0), coordRind(0),virtElmRind(0), 
  nVrtxSln(0), nCellSln(0), nVrtxFld(0), nCellFld(0)
  {
    if (cg_open(myCgFileName.c_str(), CG_MODE_WRITE, &indexFile)) cg_error_exit();
  };
  ~cgnsWriter()
  { cg_close(indexFile);};

  void setUnits(MassUnits_t mu, LengthUnits_t lu, 
                TimeUnits_t tu, TemperatureUnits_t tpu, AngleUnits_t au);
  void setBaseItrData(std::string bsitrname, int ntstp, double tval);
  void setZoneItrData(std::string zitrname, std::string grdptr, std::string slnptr);
  void setIntData(std::string intname, int intval);
  void setZone(std::string zName, ZoneType_t zt);
  void setTimestamp(std::string trimmed_base_t);
  void setNVrtx(int nVt);
  void setNCell(int nCl);
  void setGlobalNCell(int gnCl);
  void setGridXYZ(vect<double>::v1d x, vect<double>::v1d y, vect<double>::v1d z);
  void setCoordRind(int rind);  
  void setVirtElmRind(int rind); 
  void setPconnVec(const vect<int>::v1d& pConnVec); 
  void setPatchNo(int patchNo);
  void setBcflag(int bcflag);
  void setCnstrtype(int cnstr_type);
  void setVolCellFacesNumber(int nVolCellFaces);
  void setPconnGhostDescriptor(int ghostDescriptor);
  void setSection(std::string sName, ElementType_t st, vect<int>::v1d elmConn);
  void setGlobalSection(std::string sName, ElementType_t st, vect<int>::v1d elmConn);
  void setGlobalSection(std::string gsName, ElementType_t gst);
  void resetSections();
  void resetGlobalSections();
  void setSolutionNode(std::string ndeName, GridLocation_t slnLoc);
  void writeSolutionNode(std::string ndeName, GridLocation_t slnLoc, int emptyFlag, int virtFlag, int typeFlag);
  void writeSolutionField(std::string fname, std::string ndeName, DataType_t dt, void* data);
  void writeGridToFile();
  void writeZoneToFile(int volSurfFlag);
  void writeWinToFile();
  void setFluidUnitsMap();
  void setFluidDimMap();
  void setFluidMagMap();
  void setiFluidUnitsMap();
  void setiFluidDimMap();
  void setiFluidMagMap();
  //void writeVolGs(double* data);

private:
  // general CGNS manager
  std::string myCgFileName;
  ZoneType_t zoneType;
  ElementType_t sectionType, gsectionType;
  int rmin[3], rmax[3];
  cgsize_t cgCoreSize[9];
  std::string baseName, zoneName, sectionName, gsectionName;
  int indexFile, indexBase;
  int indexZone, indexCoord;
  int indexSection;
  int cellDim, physDim;
  std::string trimmed_base_t;
  int nBase;
  int nZone;
  int nVrtx;
  int nCell;
  int gnCell;
  int nSection;
  int gnSection;
  int nVrtxElem;
  int coordRind;
  int virtElmRind;
  int pConnGhostDescriptor;
  vect<int>::v1d pConnVec;
  vect<int>::v1d elmConn;
  int patchNo;
  int bcflag;
  int cnstr_type;
  int intVal;
  int nVolCellFaces;
  // units
  MassUnits_t massU;
  LengthUnits_t lengthU;
  TimeUnits_t timeU;
  TemperatureUnits_t tempU;
  AngleUnits_t angleU;
  // time-dependent data
  std::string baseItrName, zoneItrName;
  std::string intName;
  std::string gridCrdPntr, flowSlnPntr;
  int nTStep;
  double timeLabel;
  // mesh data
  vect<double>::v1d xCrd, yCrd, zCrd;
  vect<std::string>::v1d zoneNames;
  vect<ZoneType_t>::v1d zoneTypes;
  vect<std::string>::v1d sectionNames;
  vect<std::string>::v1d gsectionNames;
  vect<ElementType_t>::v1d sectionTypes;
  vect<ElementType_t>::v1d gsectionTypes;
  vect<int>::v2d elmConns;
  vect<int>::v2d gelmConns;
  vect<int>::v1d nCells;
  vect<int>::v1d gnCells;
  // solution data 
  int nSlnFld;
  int nVrtxSln;
  int nCellSln;
  int nVrtxFld;
  int nCellFld;
  std::map<std::string, int> slnNameNFld; // number of filed in the given solution name
  vect<int>::v1d solutionIdx; // solution node indices
  std::map<std::string, GridLocation_t> solutionNameLocMap;
  // maps sol node to index (eg. NodeData0.0=1,ElemData0.0=2)
  std::map<std::string, int> solutionNameSolIdxMap;  
  std::map<int,std::pair<int,keyValueList> > solutionMap; // (#sln, <slnIdx, [fldIdx, fldName]>)
  vect<std::string>::v1d solutionName;
  vect<GridLocation_t>::v1d solutionGridLocation;
  // Rocstar variable units lookup dictionary
  std::map <std::string, std::string> fluidUnitsMap;
  std::map <std::string, float*> fluidDimMap;
  std::map <std::string, int> fluidMagMap;
  std::map <std::string, std::string> ifluidUnitsMap;
  std::map <std::string, float*> ifluidDimMap;
  std::map <std::string, int> ifluidMagMap;
};

#endif
