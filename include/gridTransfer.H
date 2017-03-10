/* General purpose solution transfer class */

#ifndef __GRIDTRANSFER_H
#define __GRIDTRANSFER_H

// standard
#include <string>

// Nemosys
#include <cgnsAnalyzer.H>
#include <cgnsWriter.H>
#include <vtkAnalyzer.H>
#include <baseInterp.H>

// MAdLib headers
#include <MAdLib.h>
#include <NodalDataManager.h>

// GMsh
#include <GModel.h>

class gridTransfer: public cgnsAnalyzer {
public:
  gridTransfer(std::string srcFname, std::string trgFname);
 ~gridTransfer();

  // loading
  void loadSrcCgSeries(int nCg);
  void loadSrcCg();
  void loadTrgCg();

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
  bool zoneHasPane(cgnsAnalyzer* cgObj, int zoneIdx);
  bool paneHasPatchNo(cgnsAnalyzer* cgObj, int zoneIdx);
  int getPanePatchNo(cgnsAnalyzer* cgObj, int zoneIdx);
  int getPaneCnstrType(cgnsAnalyzer* cgObj, int zoneIdx);

  // data processing
  void stitchMe(cgnsAnalyzer* cgObj, int zoneIdx, int verb=0);
  void stitchFldBc(cgnsAnalyzer* cgObj, int zoneIdx, int verb=0); 
  void stitchMe(gridTransfer* cgObj);

  // transfer actions
  void exportMeshToMAdLib(std::string gridName);
  void exportMeshBndToMAdLib(std::string gridName);
  void convertToMsh(std::string gridName);
  void convertToSTL(std::string gridName, std::string prefix);
  void convertToVTK(std::string gridName, bool withSolution = false, std::string dist = "");
  void exportNodalDataToMAdLib();
  void exportToGModel(std::string msh);
  void transfer();
  double calcTransAcc(std::string slnName); 

  // quality reporting
  void gridStats();
  void gridCheck();
  void gridHist();

  // writing actions
  void writeTrgCg(std::string cgFName);
  
  // helper methods
  int getElmIdx(std::string msh, std::vector<double>& xyz);
  int getBaryCrds(std::string msh, std::vector<double>& xyz, 
                  std::vector<double>& baryCrds, std::vector<int>& vrtIds);

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
  bool isTransferred;
  
  // mesh data
private:
  MAd::pGModel srcModel;
  MAd::pGModel trgModel;
  MAd::pMesh srcMesh;
  MAd::pMesh trgMesh;

  // Gmsh data
public:
  GModel* srcGModel;
  GModel* trgGModel;

  // BCs data (example implementation)
private:
  std::vector<int> bcFlag;
  std::vector<int> bFlag;
  std::vector<int> patchNo; 
};

#endif

