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

    public:
      bool isCompatible() const {return _isSupported;}

      int getNumberOfNode() const {return numNdes;}
      int getNumberOfElement() const {return numElms;}
      int getNumberOfNodeSet() const {return _ndeSet.size();}
      int getNumberOfElementBlock() const {return _elmBlock.size();}
      int getNumberOfSideSets() const {return _sideSet.size();}

      void addNdeSet(ndeSetType& ns) { _ndeSet.push_back(ns); }
      void addElmBlk(elmBlockType& eb) { _elmBlock.push_back(eb); }
      
      std::string getBlockName(int id) const;
      elementType getBlockElmType(int id) const;

      void write();

   private:
      // populate quanitites needed to write
      // pntmesh file
      void exoPopulate();
      void wrnErrMsg(int errCode, std::string errMsg = "None");

    private:
      std::string _ifname;
      bool _isSupported; // false is non Tri/Tet elements were found
      bool _isPopulated;
      bool _isOpen;
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
      std::vector<char*> ndeSetNames;
      std::vector<char*> elmBlockNames;
      std::vector<char*> sideSetNames;
      std::vector<int> elmBlockOffset;
};

} // end namespace exoMesh
#endif //exoMesh
