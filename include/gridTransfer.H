/* General purpose solution transfer class */

#ifndef __GRIDTRANSFER_H
#define __GRIDTRANSFER_H

// Nemosys
#include "cgnsAnalyzer.H"
#include <string>

// MAdLib headers
#include "MAdLib.h"

class gridTransfer: public cgnsAnalyzer {
public:
  gridTransfer(std::string srcFname, std::string trgFname);
 ~gridTransfer();

  // loading
  void loadSrcCgSeries(int nCg);
  void loadTrgCgSeries(int nCg);

  // cgObjects access
  int getNCgObj();

  // name/misc access
  std::string getBaseName();
  std::string getBaseName(int indx);
  std::string getBaseNameSrc();
  std::string getBaseNameTrg();
  std::string getCgFName(int indx);
  std::string getBaseItrName(int indx);
  int getNTStep(int indx);
  double getTimeStep(int indx);
  std::string getZoneItrName(int indx, int zidx);
  std::string getGridCrdPntr(int indx, int zidx);
  std::string getSolutionPntr(int indx, int zidx);

  // zone data access
  int getNZone(int indx);
  std::string getZoneName(cgnsAnalyzer* cgObj, int zoneIdx);
  std::string getZoneName(int cgIdx, int zoneIdx);
  ZoneType_t getZoneType(int indx, int zidx);
  std::string getSectionName(int cgIdx, int zoneIdx);
  int getElementType(int cgIdx, int zoneIdx);
  int getZoneNVrtx(cgnsAnalyzer* cgObj, int zoneIdx);
  int getZoneNCell(cgnsAnalyzer* cgObj, int zoneIdx);
  std::vector<double> getZoneCoords(cgnsAnalyzer* cgObj, int zoneIdx, int dim); 
  std::vector<int> getZoneRealConn(cgnsAnalyzer* cgObj, int zoneIdx); 
  int getZoneRealSecType(cgnsAnalyzer* cgObj, int zoneIdx);

  // pane data access
  int getPaneBcflag(cgnsAnalyzer* cgObj, int zoneIdx);
  bool paneHasPatchNo(cgnsAnalyzer* cgObj, int zoneIdx);
  int getPanePatchNo(cgnsAnalyzer* cgObj, int zoneIdx);
  int getPaneCnstrType(cgnsAnalyzer* cgObj, int zoneIdx);

  // data processing
  void stitchMe(cgnsAnalyzer* cgObj, int zoneIdx);
  void stitchFldBc(cgnsAnalyzer* cgObj, int zoneIdx); 
  void stitchMe(gridTransfer* cgObj);

  // transfer actions
  void dummy();
  void exportToMAdLib(std::string gridName);
  void convertToMsh(std::string gridName);
  void convertToVtk(std::string gridName, bool withSolution = false);


  // management data
private:
  std::string srcCgFName;
  std::string trgCgFName;
  std::string baseCgFNameSrc;
  std::string baseCgFNameTrg;
  std::vector<std::string> srcCgFNames;
  std::vector<std::string> trgCgFNames;
  std::vector<std::string> nowCgFNames;
  int padSizeSrc, padSizeTrg;
  std::vector<cgnsAnalyzer*> srcCgObjs;
  std::vector<cgnsAnalyzer*> trgCgObjs;
  std::vector<cgnsAnalyzer*> nowCgObjs;
  
  // mesh data
private:
  MAd::pGModel srcModel;
  MAd::pGModel trgModel;
  MAd::pMesh srcMesh;
  MAd::pMesh trgMesh;

  // BCs data (example implementation)
private:
  std::vector<int> bcFlag;
  std::vector<int> bFlag;
  std::vector<int> patchNo; 
};

#endif

