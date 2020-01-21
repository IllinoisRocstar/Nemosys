#ifndef NEMOSYS_EXOMESH_H_
#define NEMOSYS_EXOMESH_H_

#include <string>
#include <vector>

#include <vtkCellTypes.h>

#include "nemosys_export.h"

namespace NEM {
namespace MSH {
namespace EXOMesh {

/**
 * Logging method.
 * @param errCode If <0, then error and throw. If =0, continue. If >0, warning.
 * @param errMsg Message to output explaining error.
 */
void NEMOSYS_EXPORT wrnErrMsg(int errCode, const std::string &errMsg = "None");

/**
 * EXODUS Element Types.
 */
enum NEMOSYS_EXPORT elementType { TRIANGLE, QUAD, TETRA, HEX, WEDGE, OTHER };

/**
 * Boundary condition types
 */
enum NEMOSYS_EXPORT surfaceBCTag { FIXED, SYMMX, SYMMY, SYMMZ };

/**
 * Stores node set information.
 */
struct NEMOSYS_EXPORT ndeSetType {
  int id;                   // written
  int nNde;                 // written
  std::string name;         // written
  std::vector<int> ndeIds;  // written
  int ndeIdOffset;
};

/**
 * Stores element block information.
 */
struct NEMOSYS_EXPORT elmBlkType {
  int id;                 // written
  std::string name;       // written
  elementType eTpe;       // written
  int ndePerElm;          // written
  int nElm;               // written
  std::vector<int> conn;  // written
  int ndeIdOffset;
  std::vector<int> elmIds;
};

/**
 * Stores side set information.
 */
struct NEMOSYS_EXPORT sdeSetType {
  int id;                   // written
  std::string name;         // written
  int nSde;                 // written
  std::vector<int> elmIds;  // written
  std::vector<int> sdeIds;  // written
  int elmIdOffset;
};

typedef std::vector<ndeSetType> ndeSets;
typedef std::vector<elmBlkType> elmBlks;
typedef std::vector<sdeSetType> sdeSets;

/**
 * Convert EXODUS element type to VTK cell type
 * @param et EXODUS element type
 * @return VTK cell type
 */
VTKCellType NEMOSYS_EXPORT e2vEMap(elementType et);
/**
 * Convert VTK cell type to EXODUS element type
 * @param vt VTK cell type
 * @return EXODUS element type
 */
elementType NEMOSYS_EXPORT v2eEMap(VTKCellType vt);
/**
 *
 * @param tag
 * @return
 */
// TODO: documentation
surfaceBCTag NEMOSYS_EXPORT bcTagNum(std::string &tag);
/**
 *
 * @param tag
 * @return
 */
// TODO: documentation
std::string NEMOSYS_EXPORT bcTagStr(int tag);
/**
 * Convert string to EXODUS element type
 * @param tag string tag of an element type
 * @return EXODUS element type
 */
elementType NEMOSYS_EXPORT elmTypeNum(std::string tag);
/**
 * Convert EXODUS element type to string tab
 * @param et EXODUS element type
 * @return string tag of an element type
 */
std::string NEMOSYS_EXPORT elmTypeStr(elementType et);
/**
 * Get number of nodes given EXODUS element type and order
 * @param et EXODUS element type
 * @param order element order
 * @return number of nodes
 */
int NEMOSYS_EXPORT elmNumNde(elementType et, int order);
/**
 * Get number of surfaces given EXODUS element type
 * @param et EXODUS element type
 * @return number of surfaces
 */
int NEMOSYS_EXPORT elmNumSrf(elementType et);

/**
 * @brief A complete I/O class for EXODUS II file format. 
 *
 * Methods for definition, addition, removal, and modification of the elements,
 * nodes, element blocks, node sets, and side sets are implemented. The
 * implementation relies on standard EXODUS II library provided in system 
 * repository.
 */
class NEMOSYS_EXPORT exoMesh {
 public:
  exoMesh();
  /**
   * Construct exoMesh with initial file name
   * @param ifname file name
   */
  explicit exoMesh(std::string ifname);

  exoMesh(const exoMesh &) = default;
  exoMesh(exoMesh &&) = default;

  exoMesh &operator=(const exoMesh &) = default;
  exoMesh &operator=(exoMesh &&) = default;

  ~exoMesh();

  // fast access
 public:
  // TODO: this flag is no longer supported and will be deprecated soon
  bool isCompatible() const { return _isSupported; }

  /**
   * Returns total number of nodes
   */
  int getNumberOfNodes() const { return _numNdes; }

  /**
   * Returns total number of elements
   */
  int getNumberOfElements() const { return _numElms; }

  /**
   * Returns total number of element blocks
   */
  int getNumberOfElementBlocks() const { return _elmBlks.size(); }

  /**
   * Returns total number of node sets
   */
  int getNumberOfNodeSets() const { return _ndeSets.size(); }

  /**
   * Returns total number of side sets
   */
  int getNumberOfSideSets() const { return _sdeSets.size(); }

  /**
   * Returns the name the element block for given index
   * @param idx index of block
   */
  const std::string &getElmBlkName(int idx) const { return _elmBlks[idx].name; }

  /**
   * Returns the name the nodeset for given index
   * @param idx index of block
   */
  const std::string &getNdeSetName(int idx) const { return _ndeSets[idx].name; }

  /**
   * Returns the name the sideset for given index
   * @param idx index of block
   */
  const std::string &getSdeSetName(int idx) const { return _sdeSets[idx].name; }

  /**
   * Returns the id for the element block for given index
   * @param idx index of block
   */
  int getElmBlkId(int idx) const { return _elmBlks[idx].id; }

  /**
   * Returns the id for the nodeset for given index
   * @param idx index of block
   */
  int getNdeSetId(int idx) const { return _ndeSets[idx].id; }

  /**
   * Returns the id for the sideset for given index
   * @param idx index of block
   */
  int getSdeSetId(int idx) const { return _sdeSets[idx].id; }

  /**
   * Returns the element type for the block
   * @param idx index of block
   */
  elementType getElmBlkType(int idx) const { return _elmBlks[idx].eTpe; };

  /**
   * Returns the element type for the block
   * @param ebName name of the element block
   */
  elementType getElmBlkType(std::string ebName) const;

  /**
   * Returns the names of registered element blocks
   */
  const std::vector<std::string> &getElmBlkNames() const {
    return _elmBlkNames;
  };

  /**
   * Returns the names of registered node sets
   */
  const std::vector<std::string> &getNdeSetNames() const {
    return _ndeSetNames;
  };

  /**
   * Returns the names of registered side sets
   */
  const std::vector<std::string> &getSdeSetNames() const {
    return _sdeSetNames;
  };

  /**
   * Returns problem dimension
   */
  int getDimension() const { return _numDim; }

  // fast modification
 public:

  /**
   * Sets/changes block name
   * @param idx index of block
   * @param name name of the block
   */
  void setBlockName(int idx, const std::string &name) {
    _elmBlks[idx].name = name;
  }

  /**
   * Sets/changes the problem dimensionality
   * @param dim dimension of problem 2D/3D
   */
  void setDimension(int dim) {
    _numDim = dim;
  }

  // adding nodes, element blocks, node sets, and side sets
 public:
  /**
   * Add nodes to the database
   * @param x
   * @param y
   * @param z
   */
  void addNde(double x, double y, double z) {
    _xCrds.emplace_back(x);
    _yCrds.emplace_back(y);
    _zCrds.emplace_back(z);
    ++_numNdes;
  }
  /**
   * Add nodes to the database
   * @param xyz 3-vector of form (x, y, z)
   */
  void addNde(const std::vector<double> &xyz) {
    addNde(xyz[0], xyz[1], xyz[2]);
  }
  /**
   * Add element block to the database
   * @param eb EXODUS element block
   */
  void addElmBlk(const elmBlkType &eb) {
    _elmBlks.emplace_back(eb);
    _isPopulated = false;
  }
  /**
   * Add node set to the database
   * @param ns EXODUS node set
   */
  void addNdeSet(const ndeSetType &ns) {
    _ndeSets.emplace_back(ns);
    _isPopulated = false;
  }
  /**
   * Add side set to the database
   * @param ss EXODUS side set
   */
  void addSdeSet(const sdeSetType &ss) {
    _sdeSets.emplace_back(ss);
    _isPopulated = false;
  }

  // removing nodes, element blocks, node sets, and side sets
 public:
  /**
   * Remove an element block by name
   * @note Calls exoPopulate.
   * @param blkName element block name
   */
  void removeElmBlkByName(const std::string &blkName);

  // finding element blocks
 public:
  /**
   * Finds index of the first block containing element
   * @param elmId element id
   * @return element block index
   */
  int findElmBlkIdxByElmId(int elmId) const;  // returns -1 if not found
  /**
   * Finds index of the block containing most or all of the list (slower)
   * @note returns -1 if not found
   * @param elmIds vector of element ids
   * @return element block index
   */
  int findElmBlkIdxByElmIdLst(const std::vector<int> &elmIds) const;

  // finding elements
 public:
  /**
   * Finds all elements that are within the block and generates a list of them
   * @param blkIdx element block index
   * @param elmIds vector of element ids
   * @param allIn true if all elements are in given element block
   * @return subset of given element ids in given element block
   */
  std::vector<int> lstElmInBlk(int blkIdx, const std::vector<int> &elmIds,
                               bool &allIn) const;

  // augmentation and maintenance
 public:
  /**
   * Removes a list of elements from an element block
   * @param blkIdx element block index
   * @param idLst vector of element ids
   */
  // TODO: Does NOT update side sets properly!
  void removeByElmIdLst(int blkIdx, const std::vector<int> &idLst);
  
  /**
   * Creates a new element block and augments previous owners
   * @param name new element block name
   * @param idLst vector of element ids
   */
  void addElmBlkByElmIdLst(const std::string &name, std::vector<int> &idLst);

  /**
   * Creates a new node set and augments previous ones if needed
   * @param name new node set name
   * @param idLst vector of node ids
   */
  void addNdeSetByNdeIdLst(const std::string &name,
                           const std::vector<int> &idLst);
  /**
   * Filter nodal coordinates and snap to zero
   * @param tol tolerance used to determine if zero
   */
  void snapNdeCrdsZero(double tol = 1e-5);

  /**
   * Merges duplicated and nodes within given proximity
   * @param tol tolerance used in determining if a line intersects a surface
   */
  void mergeNodes(double tol = 1e-15);

  /**
   * scales the nodal coordinates 
   * @param sc scale factor
   */
  void scaleNodes(double sc = 1.0);

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
  /**
   * Write to file name specified at construction or using setFileName method
   */
  void write();
  /**
   * Resets the class and reads from file name provided
   */
  void read(const std::string &ifname = std::string());
  /**
   * Resets the EXODUS database. Clears nodes, element blocks, node sets, side
   * sets.
   */
  void reset();

  // Miscellaneous
 public:
  /**
   * Print out a report of the current EXODUS database.
   */
  void report() const;

  /**
   * Toggle verbosity.
   */
  void togVrb() { _isVerbose = !_isVerbose; };

 private:
  // populate quantities needed to write EXODUS file
  void exoPopulate(bool updElmLst = false);

  // EXODUS data
 private:
  // Dimension
  int _numDim;

  // Nodes
  int _numNdes;
  std::vector<double> _xCrds;
  std::vector<double> _yCrds;
  std::vector<double> _zCrds;

  // Element Blocks
  int _numElms;
  elmBlks _elmBlks;
  std::vector<std::string> _elmBlkNames;

  // Node Sets
  ndeSets _ndeSets;
  std::vector<std::string> _ndeSetNames;

  // Side Sets
  sdeSets _sdeSets;
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
}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_EXOMESH_H_
