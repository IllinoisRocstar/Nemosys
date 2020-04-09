#ifndef __CGNSANALYZER_H
#define __CGNSANALYZER_H

// Nemosys headers
#include "nemosys_export.h"
#include "meshBase.H"

// VTK
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>

// MAdLib
#include <MAdLib.h>

// CGNS
#include <cgnslib.h>

// ANN
#include <ANN/ANN.h>

// Eigen
//#include <Dense>

// others
#include <vector>
#include <string>
#include <map>

// enumerations
enum solution_type_t {
  NODAL, ELEMENTAL, UNKNOWN
};

// typedefs
typedef std::map<int, std::string> keyValueList;
//typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixInt;
//typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorInt;

class NEMOSYS_EXPORT solutionData
{
//TODO change vtk specific stuff to use meshBase
    typedef std::vector<double> vecSlnType;

  public:
    solutionData(std::string sName, solution_type_t inDataType)
        : dataName(std::move(sName)), nDim(0), nData(0), dataType(inDataType) {}

    ~solutionData() = default;

    void appendData(const vecSlnType &data, int inNData, int inNDim);
    void getData(vecSlnType &inBuff, int &outNData, int &outNDim) const;
    void getData(vecSlnType &inBuff, int &outNData, int &outNDim,
                 const std::vector<bool> &mask) const;
    void rmvDataIdx(const std::vector<int> &rmvIdx);

    solutionData *getPtr() { return this; };

    std::string getDataName() const { return dataName; };

    int getNDim() const { return nDim; };

    int getNData() const { return nData; };

    solution_type_t getDataType() const { return dataType; };

    void clearBuff()
    {
      dataName = "";
      dataType = UNKNOWN;
      nDim = 0;
      nData = 0;
      slnData.clear();
    };

  private:
    vecSlnType slnData;
    std::string dataName;
    int nDim;
    int nData;
    solution_type_t dataType;
};


class NEMOSYS_EXPORT cgnsAnalyzer
{
  public:
    explicit cgnsAnalyzer(std::string fname, int verb = 0) :
        cgFileName(std::move(fname)), isUnstructured(false),
        zoneType(CGNS_ENUMV(ZoneTypeNull)),
        indexFile(-1), indexBase(-1), indexZone(-1), indexCoord(-1),
        cellDim(0), physDim(0), nVertex(0), nRindNdeStr(0), nElem(0),
        nVrtxElem(0),
        solutionDataPopulated(false),
        kdTree(nullptr), kdTreeElem(nullptr), vrtxCrd(nullptr), vrtxIdx(nullptr),
        searchEps(1e-9),
        isMltZone(false), vtkMesh(nullptr), _verb(verb), _rindOff(false)
    {
      cgRindCellIds.clear();
      cgRindNodeIds.clear();
    };

    // destructor
    virtual ~cgnsAnalyzer()
    {
      // close file if still open
      if (indexFile > 0)
        cg_close(indexFile);
      delete kdTree;
      delete kdTreeElem;
      if (vrtxCrd)
        annDeallocPts(vrtxCrd);
      if (vrtxIdx)
        annDeallocPts(vrtxIdx);
      clearAllSolutionData();
    };

    // top-level loading mesh
    void loadGrid(const std::string &fname, int verb = 0)
    {
      cgFileName = fname;
      loadGrid(verb);
    };

    // closing cgns file
    virtual void closeCG()
    {
      if (indexFile > 0)
      {
        cg_close(indexFile);
        indexFile = -1;
      }
    };

    // general purpose IO
    void loadGrid(int verb = 0);
    void loadZone(int zIdx, int verb = 0);

    // mesh information access
    int getIndexFile();
    int getIndexBase();
    int getCellDim();
    std::string getFileName();
    std::string getBaseName();
    std::string getZoneName();
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
    CGNS_ENUMT(MassUnits_t) getMassUnit();
    CGNS_ENUMT(LengthUnits_t) getLengthUnit();
    CGNS_ENUMT(TimeUnits_t) getTimeUnit();
    CGNS_ENUMT(TemperatureUnits_t) getTemperatureUnit();
    CGNS_ENUMT(AngleUnits_t) getAngleUnit();
    CGNS_ENUMT(ZoneType_t) getZoneType();
    bool isStructured();
    std::vector<double> getVertexCoords();
    std::vector<double> getVertexCoords(int vrtxId);
    double getVrtXCrd(int vrtxId);
    std::vector<double> getVrtXCrd();
    double getVrtYCrd(int vrtxId);
    std::vector<double> getVrtYCrd();
    double getVrtZCrd(int vrtxId);
    std::vector<double> getVrtZCrd();
    std::vector<cgsize_t> getElementConnectivity(int elemId);

    // other general purpose mesh query (such as virtual meshes)
    void getSectionNames(std::vector<std::string> &secNames);
    void getSectionConn(std::string secName, std::vector<cgsize_t> &conn, int &nElm);
    vtkSmartPointer<vtkDataSet> getSectionMesh(std::string secName);

    // get solution information
    void clearAllSolutionData();
    void getSolutionDataNames(std::vector<std::string> &list);
    solution_type_t
    getSolutionData(std::string sName, std::vector<double> &slnData);
    solutionData *getSolutionDataObj(std::string sName); // reads from CGNS file
    int getNVertexSolution();
    int getNCellSolution();
    solution_type_t
    getSolutionDataStitched(std::string sName, std::vector<double> &slnData,
                            int &outNData, int &outNDim);

    // append and query new solution/bc fields
    void appendSolutionData(std::string sName, std::vector<double> &slnData,
                            solution_type_t dt, int inNData, int inNDim);
    void appendSolutionData(std::string sName, double slnData,
                            solution_type_t dt, int inNData, int inNDim);
    bool delAppSlnData(std::string sName);
    void getAppendedSolutionDataName(std::vector<std::string> &appSName);

    // access solution data related management data structures
    std::vector<std::string> getSolutionNodeNames();
    std::vector<CGNS_ENUMT(GridLocation_t)> getSolutionGridLocations();
    std::map<int, std::pair<int, keyValueList> > getSolutionMap();
    std::map<std::string, CGNS_ENUMT(GridLocation_t)> getSolutionNameLocMap();

    // mesh processing
    void exportToMAdMesh(MAd::pMesh MAdMesh);
    virtual void stitchMesh(cgnsAnalyzer *inCg, bool withFields = false);
    void classifyMAdMeshOpt(MAd::pMesh MAdMesh);
    void classifyMAdMeshBnd(MAd::pMesh MAdMesh);
    void unclassifyMAdMeshBnd(MAd::pMesh MAdMesh);
    void checkVertex();
    bool checkElmConn(int nSharedNde);
    std::vector<double> getElmCntCoords(MAd::pMesh msh);

    // export to vtk format without going through madlib/gmsh (needed for simmetrix interface)
    void exportToVTKMesh();
    vtkSmartPointer<vtkDataSet> getVTKMesh();

    // write data from existing vtkdataset into open cgns file (for replacing data in file)
    void overwriteSolData(meshBase *mbObj);
    void overwriteSolData(const std::string &fname, const std::string &ndeName,
                          int slnIdx, CGNS_ENUMT(DataType_t) dt, void *data);

    // utilities
    void writeSampleStructured();
    void writeSampleUnstructured();

    // others (special purpose)
    bool isMultiZone();
    bool isCgRindNode(int cgNdeId);
    bool isCgRindCell(int cgCellId);
    void cleanRind();


  protected:
    void populateSolutionDataNames();
    void buildVertexKDTree();
    void buildElementKDTree();
    void loadSolutionDataContainer(int verb = 0);
    virtual void stitchFields(cgnsAnalyzer *inCg);
    CGNS_ENUMT(ElementType_t) getSectionType(std::string secName);

  protected:
    std::string cgFileName;
    // type of the grid
    std::string baseName, zoneName, sectionName;
    bool isUnstructured;
    CGNS_ENUMT(ZoneType_t) zoneType;
    CGNS_ENUMT(ElementType_t) sectionType;
    // units
    CGNS_ENUMT(MassUnits_t) massU;
    CGNS_ENUMT(LengthUnits_t) lengthU;
    CGNS_ENUMT(TimeUnits_t) timeU;
    CGNS_ENUMT(TemperatureUnits_t) tempU;
    CGNS_ENUMT(AngleUnits_t) angleU;
    // grid ints
    cgsize_t rmin[3], rmax[3];
    cgsize_t cgCoreSize[9];
    int indexFile, indexBase;
    int indexZone, indexCoord;
    int indexSection;
    int cellDim, physDim;
    int nBase;
    int nZone;
    int nVertex;
    int nRindNdeStr;
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
    std::vector<cgsize_t> elemConn;
    // solution data
    bool solutionDataPopulated;
    int nSolution;
    int nField;
    std::vector<solutionData *> slnDataCont;
    std::map<std::string, CGNS_ENUMT(GridLocation_t)> solutionNameLocMap;
    std::map<int, std::pair<int, keyValueList> > solutionMap; // (#sln, <slnIdx, (fldIdx, fldName)>)
    std::vector<std::string> solutionName;
    std::vector<CGNS_ENUMT(GridLocation_t)> solutionGridLocation;
    std::vector<std::string> appendedSolutionName;
    // export variables
    std::map<int, int> MAdToCgnsIds;
    std::map<int, int> cgnsToMAdIds;
    // search support data structures
    ANNkd_tree *kdTree;
    ANNkd_tree *kdTreeElem;
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
    std::vector<int> cgRindCellIds;
    std::vector<int> cgRindNodeIds;
    // other flags
    int _verb;
    bool _rindOff;
};

#endif
