#ifndef NEMOSYS_EXOMESH_H_
#define NEMOSYS_EXOMESH_H_

#include <string>
#include <vector>

#include <vtkCellTypes.h>

#include "nemosys_export.h"

namespace EXOMesh {

void NEMOSYS_EXPORT wrnErrMsg(int errCode, const std::string &errMsg = "None");

enum NEMOSYS_EXPORT elementType { TRIANGLE, QUAD, TETRA, HEX, WEDGE, OTHER };

enum NEMOSYS_EXPORT surfaceBCTag { FIXED, SYMMX, SYMMY, SYMMZ };

struct NEMOSYS_EXPORT ndeSetType {
  int id;                   // written
  int nNde;                 // written
  std::string name;         // written
  std::vector<int> ndeIds;  // written
  int ndeIdOffset;
};

struct NEMOSYS_EXPORT elmBlockType {
  int id;                 // written
  std::string name;       // written
  elementType eTpe;       // written
  int ndePerElm;          // written
  int nElm;               // written
  std::vector<int> conn;  // written
  int ndeIdOffset;
  std::vector<int> elmIds;
};

struct NEMOSYS_EXPORT sideSetType {
  int id;
  std::vector<int> elmIds;
  std::vector<int> sideIds;
};

typedef std::vector<ndeSetType> ndeSet;
typedef std::vector<elmBlockType> elmBlock;
typedef std::vector<sideSetType> sideSet;

VTKCellType NEMOSYS_EXPORT e2vEMap(elementType et);
elementType NEMOSYS_EXPORT v2eEMap(VTKCellType vt);
surfaceBCTag NEMOSYS_EXPORT bcTagNum(std::string &tag);
std::string NEMOSYS_EXPORT bcTagStr(int tag);
elementType NEMOSYS_EXPORT elmTypeNum(std::string tag);
std::string NEMOSYS_EXPORT elmTypeStr(elementType tag);
int NEMOSYS_EXPORT elmNumNde(elementType et, int order);
int NEMOSYS_EXPORT elmNumSrf(elementType et);

class NEMOSYS_EXPORT exoMesh {
 public:
  exoMesh();
  explicit exoMesh(std::string ifname);

  exoMesh(const exoMesh &) = default;
  exoMesh(exoMesh &&) = default;

  exoMesh &operator=(const exoMesh &) = default;
  exoMesh &operator=(exoMesh &&) = default;

  ~exoMesh();

  // fast access
 public:
  bool isCompatible() const { return _isSupported; }
  int getNumberOfNode() const { return _numNdes; }
  int getNumberOfElement() const { return _numElms; }
  int getNumberOfNodeSet() const { return _ndeSet.size(); }
  int getNumberOfElementBlock() const { return _elmBlock.size(); }
  int getNumberOfSideSets() const { return _sdeSet.size(); }
  const std::string &getBlockName(int id) const { return _elmBlock[id].name; }
  const std::string &getNdeSetName(int id) const { return _ndeSet[id].name; }
  int getElmBlkId(int id) const { return _elmBlock[id].id; }
  int getNdeSetId(int id) const { return _ndeSet[id].id; }
  elementType getBlockElmType(int id) const { return _elmBlock[id].eTpe; };
  void report() const;
  void togVrb() { _isVerbose = !_isVerbose; };
  const std::vector<std::string> &getNdeSetNames() const {
    return _ndeSetNames;
  };
  const std::vector<std::string> &getElmBlkNames() const {
    return _elmBlkNames;
  };

  // fast modification
 public:
  void setBlockName(int id, const std::string &name) {
    _elmBlock[id].name = name;
  }

  // adding nodes, element blocks, node sets, and side sets
 public:
  void addNde(double x, double y, double z) {
    _xCrds.emplace_back(x);
    _yCrds.emplace_back(y);
    _zCrds.emplace_back(z);
    ++_numNdes;
  }
  void addNde(const std::vector<double> &xyz) {
    addNde(xyz[0], xyz[1], xyz[2]);
  }
  void addElmBlk(const elmBlockType &eb) {
    _elmBlock.emplace_back(eb);
    _isPopulated = false;
  }
  void addNdeSet(const ndeSetType &ns) {
    _ndeSet.emplace_back(ns);
    _isPopulated = false;
  }
  void removeElmBlkByName(const std::string &blkName);

  // augmentation and maintenance
 public:
  // finds index of the first block containing element
  int findElmBlkIdx(int elmId) const;  // returns -1 if not found
  // finds index of the block containing most or all of the list (slower)
  // returns -1 if not found
  int findElmLstBlkIdx(const std::vector<int> &elmIds) const;
  // finds all elements that are within the block and generates a list of them
  // if all in, sets allIn to true
  std::vector<int> lstElmInBlck(int blkId, const std::vector<int> &elmIds,
                                bool &allIn) const;
  // removes a list of elements from an element block
  void removeByElmIdLst(int blkId, const std::vector<int> &idLst);
  // creates a new element block and augments previous owners
  void addElmBlkByElmIdLst(const std::string &name, std::vector<int> &idLst);
  // creates a new node set and augments previous ones if needed
  void addNdeSetByNdeIdLst(const std::string &name,
                           const std::vector<int> &idLst);
  // filter nodal coordinates and snap to zero
  void snapNdeCrdsZero(double tol = 1e-5);
  // merges duplicated and nodes within given proximity
  void mergeNodes(double tol = 1e-15);

 public:
  /**
   * Stitch another mesh into the current.
   *
   * The stitching appends nodes, element blocks, node sets, and side sets to
   * the end of the current exoMesh. The appended items will be re-indexed off
   * the end of the current.
   * @param otherMesh Other mesh to stitch into current mesh.
   */
  void stitch(const exoMesh &otherMesh);

  // I/O
 public:
  /**
   * Change the file name. Affects calls to the write method.
   * @param fName file name
   */
  void setFileName(const std::string &fName) { _ifname = fName; }
  void write();
  void reset();
  void read(const std::string &ifname = std::string());

 private:
  // populate quantities needed to write EXODUS file
  void exoPopulate(bool updElmLst = false);

  // EXODUS data
 private:
  // Nodes
  int _numNdes;
  std::vector<double> _xCrds;
  std::vector<double> _yCrds;
  std::vector<double> _zCrds;

  // Element Blocks
  int _numElms;
  elmBlock _elmBlock;
  std::vector<std::string> _elmBlkNames;

  // Node Sets
  ndeSet _ndeSet;
  std::vector<std::string> _ndeSetNames;

  // Side Sets
  sideSet _sdeSet;
  std::vector<std::string> _sdeSetNames;

  // Metadata
  int _fid;
  float _api_v;
  float _dbs_v;
  int _exErr;

  // exoMesh internal data
 private:
  std::string _ifname;
  bool _isSupported;  // false is non Tri/Tet elements were found
  bool _isPopulated;
  bool _isOpen;
  bool _isVerbose;
  std::vector<std::vector<int>> glbConn;
};

}  // namespace EXOMesh

#endif  // NEMOSYS_EXOMESH_H_
