#ifndef exoMesh_H
#define exoMesh_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <vtkCellTypes.h>

class meshBase;

namespace EXOMesh
{

void wrnErrMsg(int errCode, std::string errMsg = "None");

enum elementType {
  QUAD,
  TRIANGLE,
  TETRA,
  WEDGE,
  HEX,
  OTHER
};

enum surfaceBCTag {
  FIXED,
  SYMMX,
  SYMMY,
  SYMMZ
};

struct ndeSetType {
  int id;
  int nNde;
  std::string name;
  bool usrNdeIds;
  std::vector<int> ndeIds;
  std::vector<std::vector<double> > crds;
};

struct elmBlockType {
  int id;
  int nElm;
  std::string name;
  int ndePerElm;
  elementType eTpe;
  int ndeIdOffset;
  std::vector<int> elmIds;
  std::vector<int> conn;
};

struct sideSetType {
  int id;
  std::vector<int> elmIds;
  std::vector<int> sideIds;
};

typedef std::vector<ndeSetType> ndeSet;
typedef std::vector<elmBlockType> elmBlock;
typedef std::vector<sideSetType> sideSet;

VTKCellType e2vEMap (elementType et);
elementType v2eEMap (VTKCellType vt);
surfaceBCTag bcTagNum(std::string& tag);
std::string bcTagStr(int tag);
elementType elmTypeNum (std::string tag);
std::string elmTypeStr (elementType tag);
int elmNumNde (elementType et, int order);
int elmNumSrf (elementType et);

class exoMesh
{

    public:
      exoMesh(std::string ifname);
      ~exoMesh();

    // fast access
    public:
      bool isCompatible() const {return _isSupported;}
      int getNumberOfNode() const {return numNdes;}
      int getNumberOfElement() const {return numElms;}
      int getNumberOfNodeSet() const {return _ndeSet.size();}
      int getNumberOfElementBlock() const {return _elmBlock.size();}
      int getNumberOfSideSets() const {return _sideSet.size();}
      std::string getBlockName(int id) const { return(_elmBlock[id].name); }
      std::string getNdeSetName(int id) const { return(_ndeSet[id].name); }
      int getElmBlkId(int id) const { return(_elmBlock[id].id); }
      int getNdeSetId(int id) const { return(_ndeSet[id].id); }
      elementType getBlockElmType(int id) const { return(_elmBlock[id].eTpe); };
      void report() const;
      void togVrb() {_isVerbose=_isVerbose?false:true; };
      std::vector<std::string> getNdeSetNames() const { return ndeSetNames;};
      std::vector<std::string> getElmBlkNames() const { return elmBlockNames;};

    // adding sets and blocks
    public:
      void addNdeSet(ndeSetType& ns) { _ndeSet.push_back(ns); _isPopulated=false; }
      void addElmBlk(elmBlockType& eb) { _elmBlock.push_back(eb);  _isPopulated=false;}
      void removeElmBlkByName(std::string blkName);

    // augmentation and maintenance
    public:
      // finds index of the first block containing element
      int findElmBlkIdx(int elmId) const; // returns -1 if not found
      // finds index of the block containing most or all of the list (slower)
      int findElmLstBlkIdx(std::vector<int> elmIds) const; // returns -1 if not found
      // finds all elements that are within the block and generates a list of them
      // if all in, sets allIn to true
      std::vector<int> lstElmInBlck(int blkId, std::vector<int> elmIds, bool allIn = false ) const;
      // removes a list of elements from an element block
      void removeByElmIdLst(int blkId, std::vector<int>& idLst);
      // creates a new element block and augments previous oweners
      void addElmBlkByElmIdLst(std::string name, std::vector<int> idLst); 
      // creates a new node set and augments previous ones if needed
      void addNdeSetByNdeIdLst(std::string name, std::vector<int> idLst); 
      // filter nodal coordinates and snap to zero
      void snapNdeCrdsZero(double tol=1e-5);
    
     
    // I/O
    public:
      void write();
      void reset();
      void read(std::string ifname=std::string());

    private:
      // populate quanitites needed to write
      // pntmesh file
      void exoPopulate(bool updElmLst=false);

    private:
      std::string _ifname;
      bool _isSupported; // false is non Tri/Tet elements were found
      bool _isPopulated;
      bool _isOpen;
      bool _isVerbose;
      int _exErr;
      ndeSet _ndeSet;
      elmBlock _elmBlock;
      sideSet _sideSet;
      int fid;
      int numNdes; 
      int numElms; 
      std::vector<double> xCrds;
      std::vector<double> yCrds;
      std::vector<double> zCrds;
      std::vector<std::vector<int>> glbConn;
      std::vector<std::string> ndeSetNames;
      std::vector<std::string> elmBlockNames;
      std::vector<std::string> sideSetNames;
      float _api_v;
      float _dbs_v;

};

} // end namespace exoMesh
#endif //exoMesh
