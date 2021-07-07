#ifndef NEMOSYS_ROCSTARCGNS_H_
#define NEMOSYS_ROCSTARCGNS_H_

// Nemosys headers
#include "nemosys_export.h"
#include "IO/cgnsAnalyzer.H"

/* Special purpose class for Rocstar CGNS files */
class NEMOSYS_EXPORT rocstarCgns : public cgnsAnalyzer
{
  public:
    rocstarCgns(std::string fname);
    rocstarCgns(const std::vector<std::string> &fnames);
    ~rocstarCgns();

  public:
    // loading/unloading
    void loadCgSeries();
    void loadCgSeries(int nCg);
    virtual void closeCG();

    // cgObjects access
    int getNCgObj();

    // name/misc access
    std::string getBaseName();
    std::string getBaseName(int indx);
    std::string getCgFName(int indx);
    std::string getBaseItrName(int indx);
    int getNTStep(int indx);
    double getTimeStep(int indx);
    std::string getZoneItrName(int indx, int zidx);
    std::string getGridCrdPntr(int indx, int zidx);
    std::string getSolutionPntr(int indx, int zidx);

    // zone data access
    int getNZone(int indx);
    std::string getZoneName(cgnsAnalyzer *cgObj, int zoneIdx);
    std::string getZoneName(int cgIdx, int zoneIdx);
    CGNS_ENUMT(ZoneType_t) getZoneType(int indx, int zidx);
    std::string getSectionName(int cgIdx, int zoneIdx);
    int getElementType(int cgIdx, int zoneIdx);
    int getZoneNVrtx(cgnsAnalyzer *cgObj, int zoneIdx);
    int getZoneNCell(cgnsAnalyzer *cgObj, int zoneIdx);
    std::vector<double>
    getZoneCoords(cgnsAnalyzer *cgObj, int zoneIdx, int dim);
    std::vector<cgsize_t> getZoneRealConn(cgnsAnalyzer *cgObj, int zoneIdx);
    int getZoneRealSecType(cgnsAnalyzer *cgObj, int zoneIdx);

    // pane data access
    int getPaneBcflag(cgnsAnalyzer *cgObj, int zoneIdx);
    int getPanePatchNo(cgnsAnalyzer *cgObj, int zoneIdx);
    int getPaneCnstrType(cgnsAnalyzer *cgObj, int zoneIdx);

    // data processing
    void stitchGroup();
    void stitchMe(cgnsAnalyzer *cgObj, int zoneIdx);
    void stitchFldBc(cgnsAnalyzer *cgObj, int zoneIdx);
    void stitchMe(rocstarCgns *cgObj);

    void setBurnBool(bool b) { _burn = b; };

  // management data
  private:
    std::string myCgFName;
    std::string baseCgFName;
    bool _burn;
    std::vector<std::string> cgFNames;
    int padSize;
    std::vector<cgnsAnalyzer *> myCgObjs;

  // mesh data
  private:
    std::vector<int> bcFlag;
    std::vector<int> bFlag;
    std::vector<int> patchNo;
};

#endif  // NEMOSYS_ROCSTARCGNS_H_
