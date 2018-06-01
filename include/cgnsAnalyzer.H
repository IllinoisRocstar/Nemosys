#ifndef __CGNSANALYZER_H
#define __CGNSANALYZER_H

// others
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <sys/time.h>
#include <math.h>

//vtk
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>

//MAd Lib
#include <ModelInterface.h>
#include <MeshDataBaseInterface.h>

// third party
#include <cgnslib.h>
#include <ANN.h>
//#include <Dense>

//#if CGNS_VERSION < 3100
//# define cgsize_t int
//#else
//# if CG_BUILD_SCOPE
//#  error enumeration scoping needs to be off
//# endif
//#endif

// enumerations
enum solution_type_t {NODAL, ELEMENTAL, UNKNOWN};

// typedefs
typedef std::map<int,std::string> keyValueList;
//typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixInt;
//typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorInt;

class solutionData {

typedef std::vector<double> vecSlnType;

public:
  solutionData(std::string sName, solution_type_t inDataType) :
  dataName(sName), dataType(inDataType), nData(0), nDim(0)
  {};

  ~solutionData()
  {};

  void appendData(const vecSlnType& data, int inNData, int inNDim);
  void getData(vecSlnType& inBuff, int& outNData, int& outNDim);
  void getData(vecSlnType& inBuff, int& outNData, int& outNDim, std::vector<bool> mask);
  solutionData* getPtr() { return(this); };
  std::string getDataName() { return(dataName); };
  int getNDim() { return(nDim); };
  int getNData() { return(nData); };
  solution_type_t getDataType() { return(dataType); };
  void clearBuff() {dataName = ""; dataType = UNKNOWN; nDim=0; nData=0; slnData.clear();};

private:
  vecSlnType slnData;
  std::string dataName;
  int nDim;
  int nData;
  solution_type_t dataType;
};

class cgnsAnalyzer {
public:
  cgnsAnalyzer(std::string fname) :
    cgFileName(fname), isUnstructured(false), zoneType(ZoneTypeNull),
    indexFile(-1), indexBase(-1), indexZone(-1), indexCoord(-1),
    nVertex(0), nElem(0), cellDim(0), physDim(0), nVrtxElem(0),
    solutionDataPopulated(false),
    searchEps(1e-9), kdTree(NULL), kdTreeElem(NULL),vrtxCrd(NULL),vrtxIdx(NULL),
    isMltZone(false),vtkMesh(0)
  {};

  ~cgnsAnalyzer()
  {
    // close file if still open
    if (indexFile>0)
       cg_close(indexFile);
    if (kdTree)
      delete kdTree;
    if (kdTreeElem)
      delete kdTreeElem;
    if (vrtxCrd)
      annDeallocPts(vrtxCrd);
    if (vrtxIdx)
      annDeallocPts(vrtxIdx);
    clearAllSolutionData(); 
  };

  // loading mesh
  void loadGrid(std::string fname, int verb=0)
  {
   cgFileName = fname;
   loadGrid(verb);
  };
  void loadGrid(int verb=0);
  void loadZone(int zIdx, int verb=0);
 
  // mesh information access
  int getIndexFile();
  int getIndexBase();
  int getCellDim();
  std::string getFileName();
  std::string getBaseName();
  virtual std::string getZoneName();
  std::string getZoneName(int nCgFile);
  std::string getSectionName();
  std::string getBaseItrName();
  int getNZone();
  int getNTStep();
  double getTimeStep();
  std::string getZoneItrName();
  std::string getGridCrdPntr();
  std::string getSolutionPntr();
  int getNVertex();
  int getNElement();
  int getPhysDim();
  int getElementType();
  int getNVrtxElem();
  MassUnits_t getMassUnit();
  LengthUnits_t getLengthUnit();
  TimeUnits_t getTimeUnit();
  TemperatureUnits_t getTemperatureUnit();
  AngleUnits_t getAngleUnit();
  ZoneType_t getZoneType();
  bool isStructured();
  std::vector<double> getVertexCoords();
  std::vector<double> getVertexCoords(int vrtxId);
  double getVrtXCrd(int vrtxId);
  std::vector<double> getVrtXCrd();
  double getVrtYCrd(int vrtxId);
  std::vector<double> getVrtYCrd();
  double getVrtZCrd(int vrtxId);
  std::vector<double> getVrtZCrd();
  std::vector<int> getElementConnectivity(int elemId);

  // get solution information
  void clearAllSolutionData();
  void getSolutionDataNames(std::vector<std::string>& list );
  solution_type_t getSolutionData(std::string sName, std::vector<double>& slnData);
  solutionData* getSolutionDataObj(std::string sName); // reads from CGNS file
  int getNVertexSolution();
  int getNCellSolution();
  solution_type_t getSolutionDataStitched(std::string sName, std::vector<double>& slnData, 
                                          int& outNData, int& outNDim );

  // append and query new solution/bc fields
  void appendSolutionData(std::string sName, std::vector<double>& slnData, 
                               solution_type_t dt, int inNData, int inNDim );
  void appendSolutionData(std::string sName, double slnData, 
                               solution_type_t dt, int inNData, int inNDim );
  bool delAppSlnData(std::string sName);
  void getAppendedSolutionDataName(std::vector<std::string>& appSName);

  // access solution data related management data structures
  std::vector<std::string> getSolutionNodeNames();
  std::vector<GridLocation_t> getSolutionGridLocations();
  std::map<int,std::pair<int,keyValueList> > getSolutionMap();
  std::map<std::string, GridLocation_t> getSolutionNameLocMap();

  // mesh processing
  void exportToMAdMesh(MAd::pMesh MAdMesh);
  virtual void stitchMesh(cgnsAnalyzer* inCg, bool withFields = false);
  void classifyMAdMeshOpt(const MAd::pMesh MAdMesh);
  void classifyMAdMeshBnd(const MAd::pMesh MAdMesh);
  void unclassifyMAdMeshBnd(const MAd::pMesh MAdMesh);
  void checkVertex();
  bool checkElmConn(int nSharedNde);
  std::vector<double> getElmCntCoords(MAd::pMesh msh);
  // export to vtk format without going through madlib/gmsh (needed for simmetrix interface)
  void exportToVTKMesh(); 
  vtkSmartPointer<vtkDataSet> getVTKMesh();

  // utilities
  void writeSampleStructured();
  void writeSampleUnstructured();

  // others (special purpose)
  bool isMultiZone();

protected:
   void populateSolutionDataNames();
   void buildVertexKDTree();
   void buildElementKDTree();
   void loadSolutionDataContainer(int verb = 0);   
   virtual void stitchFields(cgnsAnalyzer* inCg);


protected:
  std::string cgFileName;
  // type of the grid
  std::string baseName, zoneName, sectionName;
  bool isUnstructured;
  ZoneType_t zoneType;
  ElementType_t sectionType;
  // units
  MassUnits_t massU;
  LengthUnits_t lengthU;
  TimeUnits_t timeU;
  TemperatureUnits_t tempU;
  AngleUnits_t angleU;
  // grid ints
  int rmin[3], rmax[3];
  int cgCoreSize[9];
  int indexFile, indexBase;
  int indexZone, indexCoord;
  int indexSection;
  int cellDim, physDim;
  int nBase;
  int nZone;
  int nVertex;
  int nElem;
  int nSection;
  int nVrtxElem;
  // time-dependent data
  std::string baseItrName, zoneItrName;
  std::string gridCrdPntr, flowSlnPntr;
  int nTStep;
  double timeLabel;
  // mesh data
  std::vector<double> xCrd, yCrd, zCrd;
  std::vector<int> elemConn;
  // solution data
  bool solutionDataPopulated;
  int nSolution;
  int nField;
  std::vector<solutionData*> slnDataCont;
  std::map<std::string, GridLocation_t> solutionNameLocMap;
  std::map<int,std::pair<int,keyValueList> > solutionMap; // (#sln, <slnIdx, (fldIdx, fldName)>)
  std::vector<std::string> solutionName;
  std::vector<GridLocation_t> solutionGridLocation;
  std::vector<std::string> appendedSolutionName;
  // export variables
  std::map<int,int> MAdToCgnsIds;
  std::map<int,int> cgnsToMAdIds;
  // search support data structures
  ANNkd_tree* kdTree;
  ANNkd_tree* kdTreeElem;
  ANNpointArray vrtxCrd;
  ANNpointArray vrtxIdx;
  double searchEps;
  // stitching data support
  std::vector<std::string> zoneNames;
  std::vector<bool> vrtDataMask;
  std::vector<bool> elmDataMask;
  // special purpose cgns files
  bool isMltZone;
  std::vector<std::string> multZoneNames;
  vtkSmartPointer<vtkDataSet> vtkMesh;

};

#endif
