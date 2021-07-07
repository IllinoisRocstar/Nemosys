#ifndef NEMOSYS_CGNSWRITER_H_
#define NEMOSYS_CGNSWRITER_H_

// Nemosys headers
#include "nemosys_export.h"
#include "IO/cgnsAnalyzer.H"  // for keyValueList typedef

// third party
#include <cgnslib.h>

// others
#include <string>
#include <vector>
#include <map>

// typedefs
template<class T>
struct vect {
  typedef std::vector<T> v1d;
  typedef std::vector<std::vector<T>> v2d;
};

/* a class for writing CGNS files */
// TODO: This writer is heavily specific to Rocstar. This is not the original
//       intent for this class and it should be generalized to base class and 
//       a derived one for Rocstar.
class NEMOSYS_EXPORT cgnsWriter {
  public:
    cgnsWriter(std::string fname, std::string baseName, int cellDim,
               int physDim, int mode) :
        myCgFileName(std::move(fname)), baseName(std::move(baseName)),
        cellDim(cellDim), physDim(physDim), nBase(1), nZone(0),
        nVrtx(0), nCell(0),  nSection(0), gnSection(0),
        coordRind(0), virtElmRind(0),
        nSlnFld(0), nVrtxSln(0), nCellSln(0), nVrtxFld(0), nCellFld(0)
    {
      if (mode == 0) {
        if (cg_open(myCgFileName.c_str(), CG_MODE_WRITE,
                    &indexFile))
          cg_error_exit();
      } else if (mode == 1) {
        if (cg_open(myCgFileName.c_str(), CG_MODE_MODIFY,
                    &indexFile))
          cg_error_exit();
      }
      int file_type;
      cg_get_file_type(indexFile, &file_type);
    };

    ~cgnsWriter() { cg_close(indexFile); };

    void deleteFile();
    void setUnits(CGNS_ENUMT(MassUnits_t) mu, CGNS_ENUMT(LengthUnits_t) lu,
                  CGNS_ENUMT(TimeUnits_t) tu, CGNS_ENUMT(TemperatureUnits_t) tpu,
                  CGNS_ENUMT(AngleUnits_t) au);
    void setBaseItrData(std::string bsitrname, int ntstp, double tval);
    void setZoneItrData(std::string zitrname, std::string grdptr,
                        std::string slnptr);
    void setIntData(std::string intname, int intval);
    void setZone(std::string zName, CGNS_ENUMT(ZoneType_t) zt);
    int getNSections();
    void setTimestamp(std::string trimmed_base_t);
    void setNVrtx(int nVt);
    void setNCell(int nCl);
    void setGlobalNCell(int gnCl);
    void
    setGridXYZ(vect<double>::v1d x, vect<double>::v1d y, vect<double>::v1d z);
    void setCoordRind(int rind);
    void setVirtElmRind(int rind);
    void setPconnVec(const vect<int>::v1d &pConnVec);
    void setPconnLimits(int pconnProcMin, int pconnProcMax);
    void setPatchNo(int patchNo);
    void setBcflag(int bcflag);
    void setCnstrtype(int cnstr_type);
    void setVolCellFacesNumber(int nVolCellFaces);
    void setPconnGhostDescriptor(int ghostDescriptor);
    void
    setSection(std::string sName, CGNS_ENUMT(ElementType_t) st, vect<cgsize_t>::v1d elmConn);
    void setGlobalSection(std::string sName, CGNS_ENUMT(ElementType_t) st,
                          vect<int>::v1d elmConn);
    void setGlobalSection(std::string gsName, CGNS_ENUMT(ElementType_t) gst);
    void resetSections();
    void resetGlobalSections();
    void setSolutionNode(std::string ndeName, CGNS_ENUMT(GridLocation_t) slnLoc);
    void writeSolutionNode(std::string ndeName, CGNS_ENUMT(GridLocation_t) slnLoc,
                           int emptyFlag, int virtFlag);
    void
    writeSolutionField(std::string fname, std::string ndeName, CGNS_ENUMT(DataType_t) dt,
                       void *data);
    void writeGridToFile();
    void writeZoneToFile();
    void writeWinToFile();
    void setFluidUnitsMap();
    void setFluidDimMap();
    void setFluidMagMap();
    void setiFluidUnitsMap();
    void setiFluidDimMap();
    void setiFluidMagMap();
    void setBurnUnitsMap();
    void setBurnDimMap();
    void setBurnMagMap();
    void setTypeFlag(int typeFlag);
    void writeVolCellFacesNumber();

  private:
    // general CGNS manager
    std::string myCgFileName;
    CGNS_ENUMT(ZoneType_t) zoneType;
    CGNS_ENUMT(ElementType_t) sectionType, gsectionType;
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
    int typeFlag;
    vect<int>::v1d pConnVec;
    int pConnMin;
    int pConnMax;
    vect<int>::v1d elmConn;
    int patchNo;
    int bcflag;
    int cnstr_type;
    int intVal;
    int nVolCellFaces;
    // units
    CGNS_ENUMT(MassUnits_t) massU;
    CGNS_ENUMT(LengthUnits_t) lengthU;
    CGNS_ENUMT(TimeUnits_t) timeU;
    CGNS_ENUMT(TemperatureUnits_t) tempU;
    CGNS_ENUMT(AngleUnits_t) angleU;
    // time-dependent data
    std::string baseItrName, zoneItrName;
    std::string intName;
    std::string gridCrdPntr, flowSlnPntr;
    int nTStep;
    double timeLabel;
    // mesh data
    vect<double>::v1d xCrd, yCrd, zCrd;
    vect<std::string>::v1d zoneNames;
    vect<CGNS_ENUMT(ZoneType_t)>::v1d zoneTypes;
    vect<std::string>::v1d sectionNames;
    vect<std::string>::v1d gsectionNames;
    vect<CGNS_ENUMT(ElementType_t)>::v1d sectionTypes;
    vect<CGNS_ENUMT(ElementType_t)>::v1d gsectionTypes;
    vect<cgsize_t>::v2d elmConns;
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
    std::map<std::string, CGNS_ENUMT(GridLocation_t)> solutionNameLocMap;
    // maps sol node to index (eg. NodeData0.0=1,ElemData0.0=2)
    std::map<std::string, int> solutionNameSolIdxMap;
    std::map<int, std::pair<int, keyValueList> > solutionMap; // (#sln, <slnIdx, [fldIdx, fldName]>)
    vect<std::string>::v1d solutionName;
    vect<CGNS_ENUMT(GridLocation_t)>::v1d solutionGridLocation;
    // Rocstar variable units lookup dictionary
    std::map<std::string, std::string> fluidUnitsMap;
    std::map<std::string, std::vector<float>> fluidDimMap;
    std::map<std::string, int> fluidMagMap;
    std::map<std::string, std::string> ifluidUnitsMap;
    std::map<std::string, std::vector<float>> ifluidDimMap;
    std::map<std::string, int> ifluidMagMap;
    std::map<std::string, std::string> burnUnitsMap;
    std::map<std::string, std::vector<float>> burnDimMap;
    std::map<std::string, int> burnMagMap;
};

#endif  // NEMOSYS_CGNSWRITER_H_
