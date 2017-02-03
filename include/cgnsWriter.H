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
#include "GmshEntities.h"

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
  nVrtx(0), nCell(0), nSlnFld(0), 
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
  void setZone(std::string zName, ZoneType_t zt);
  void setNVrtx(int nVt);
  void setNCell(int nCl);
  void setGridXYZ(vect<double>::v1d x, vect<double>::v1d y, vect<double>::v1d z);  
  void setSection(std::string sName, ElementType_t st, vect<int>::v1d elmConn);
  void setSolutionNode(std::string ndeName, GridLocation_t slnLoc);
  void writeSolutionField(std::string fname, std::string ndeName, DataType_t dt, void* data);
  void writeGridToFile();
  void writeZoneToFile();

private:
  // general CGNS manager
  std::string myCgFileName;
  ZoneType_t zoneType;
  ElementType_t sectionType;
  int rmin[3], rmax[3];
  cgsize_t cgCoreSize[9];
  std::string baseName, zoneName, sectionName;
  int indexFile, indexBase;
  int indexZone, indexCoord;
  int indexSection;
  int cellDim, physDim;
  int nBase;
  int nZone;
  int nVrtx;
  int nCell;
  int nSection;
  int nVrtxElem;
  vect<int>::v1d elmConn;
  // units
  MassUnits_t massU;
  LengthUnits_t lengthU;
  TimeUnits_t timeU;
  TemperatureUnits_t tempU;
  AngleUnits_t angleU;
  // time-dependent data
  std::string baseItrName, zoneItrName;
  std::string gridCrdPntr, flowSlnPntr;
  int nTStep;
  double timeLabel;
  // mesh data
  vect<double>::v1d xCrd, yCrd, zCrd;
  vect<std::string>::v1d zoneNames;
  vect<ZoneType_t>::v1d zoneTypes;
  vect<std::string>::v1d sectionNames;
  vect<ElementType_t>::v1d sectionTypes;
  vect<int>::v2d elmConns;
  vect<int>::v1d nCells;
  // solution data 
  int nSlnFld;
  int nVrtxSln;
  int nCellSln;
  int nVrtxFld;
  int nCellFld;
  std::map<std::string, int> slnNameNFld; // number of filed in the given solution name
  vect<int>::v1d solutionIdx; // solution node indices
  std::map<std::string, GridLocation_t> solutionNameLocMap;
  std::map<int,std::pair<int,keyValueList> > solutionMap; // (#sln, <slnIdx, [fldIdx, fldName]>)
  vect<std::string>::v1d solutionName;
  vect<GridLocation_t>::v1d solutionGridLocation;
};

#endif
