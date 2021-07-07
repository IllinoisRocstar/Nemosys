#include "Mesh/exoMesh.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <utility>

#include <ANN/ANN.h>
#include <exodusII.h>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace MSH {
namespace EXOMesh {

VTKCellType e2vEMap(elementType et) {
  std::map<elementType, VTKCellType> eMap = {
      {elementType::TRIANGLE, VTK_TRIANGLE},
      {elementType::QUAD, VTK_QUAD},
      {elementType::TETRA, VTK_TETRA},
      {elementType::HEX, VTK_HEXAHEDRON},
      {elementType::WEDGE, VTK_WEDGE},
      {elementType::OTHER, VTK_EMPTY_CELL}};
  return eMap[et];
}

elementType v2eEMap(VTKCellType vt) {
  std::map<VTKCellType, elementType> eMap = {
      {VTK_TRIANGLE, elementType::TRIANGLE},
      {VTK_QUAD, elementType::QUAD},
      {VTK_TETRA, elementType::TETRA},
      {VTK_HEXAHEDRON, elementType::HEX},
      {VTK_WEDGE, elementType::WEDGE},
      {VTK_QUADRATIC_TRIANGLE, elementType::TRIANGLE},
      {VTK_QUADRATIC_QUAD, elementType::QUAD},
      {VTK_QUADRATIC_TETRA, elementType::TETRA},
      {VTK_QUADRATIC_HEXAHEDRON, elementType::HEX},
      {VTK_EMPTY_CELL, elementType::OTHER}};
  return eMap[vt];
}

surfaceBCTag bcTagNum(std::string &tag) {
  nemAux::toLower(tag);
  if (tag == "fixed") return surfaceBCTag::FIXED;
  if (tag == "symmx") return surfaceBCTag::SYMMX;
  if (tag == "symmy") return surfaceBCTag::SYMMY;
  if (tag == "symmz") return surfaceBCTag::SYMMZ;
  std::cerr << "Unknown surface tag " << tag << std::endl;
  throw;
}

std::string bcTagStr(int tag) {
  if (tag == FIXED) return "FIXED";
  if (tag == SYMMX) return "SYMMX";
  if (tag == SYMMY) return "SYMMY";
  if (tag == SYMMZ) return "SYMMZ";
  std::cerr << "Unknown surface tag " << tag << std::endl;
  throw;
}

elementType elmTypeNum(std::string tag) {
  nemAux::toLower(tag);
  if (tag == "triangle") return elementType::TRIANGLE;
  if (tag == "tri") return elementType::TRIANGLE;
  if (tag == "tri3") return elementType::TRIANGLE;
  if (tag == "trishell3") return elementType::TRIANGLE;

  if (tag == "quadrilateral") return elementType::QUAD;
  if (tag == "quad") return elementType::QUAD;
  if (tag == "quad4") return elementType::QUAD;
  if (tag == "shell4") return elementType::QUAD;

  if (tag == "tetrahedron") return elementType::TETRA;
  if (tag == "tetrahedral") return elementType::TETRA;
  if (tag == "tetra") return elementType::TETRA;
  if (tag == "tetra4") return elementType::TETRA;

  if (tag == "hexahedron") return elementType::HEX;
  if (tag == "hexahedral") return elementType::HEX;
  if (tag == "brick") return elementType::HEX;
  if (tag == "hex") return elementType::HEX;
  if (tag == "hex8") return elementType::HEX;

  if (tag == "wedge") return elementType::WEDGE;
  if (tag == "prismatic") return elementType::WEDGE;
  if (tag == "prism") return elementType::WEDGE;

  std::cout << "Warning : Element type " << tag << " may be unsupported.\n";
  return elementType::OTHER;
}

std::string elmTypeStr(elementType et) {
  if (et == elementType::QUAD) return "QUAD";
  if (et == elementType::TRIANGLE) return "TRIANGLE";
  if (et == elementType::QUAD) return "QUAD";
  if (et == elementType::TETRA) return "TETRA";
  if (et == elementType::HEX) return "HEX";
  if (et == elementType::WEDGE) return "WEDGE";
  return "OTHER";
}

int elmNumNde(elementType tag, int order) {
  if (tag == elementType::TRIANGLE && order == 1) return 3;
  if (tag == elementType::TRIANGLE && order == 2) return 6;
  if (tag == elementType::QUAD && order == 1) return 4;
  if (tag == elementType::QUAD && order == 2) return 9;
  if (tag == elementType::TETRA && order == 1) return 4;
  if (tag == elementType::TETRA && order == 2) return 10;
  if (tag == elementType::HEX && order == 1) return 8;
  if (tag == elementType::HEX && order == 2) return 21;
  if (tag == elementType::WEDGE && order == 1) return 6;
  if (tag == elementType::WEDGE && order == 2) return 16;
  std::cerr << "Unknown element/order combination " << tag << " " << order
            << "\n";
  throw;
}

int elmNumSrf(elementType tag) {
  if (tag == elementType::TRIANGLE) return 3;
  if (tag == elementType::QUAD) return 4;
  if (tag == elementType::TETRA) return 4;
  if (tag == elementType::HEX) return 8;
  if (tag == elementType::WEDGE) return 5;
  std::cerr << "Unknown element type " << tag << "\n";
  throw;
}

//////////////////////////////////
// exoMesh class
//////////////////////////////////

exoMesh::exoMesh()
    : _numDim(2),
      _numNdes(0),
      _numElms(0),
      _fid(-1),
      _api_v(0.0),
      _dbs_v(0.0),
      _exErr(0),
      _isSupported(true),
      _isPopulated(false),
      _isOpen(false),
      _isVerbose(false) {}

exoMesh::exoMesh(std::string ifname)
    : _numDim(2),
      _numNdes(0),
      _numElms(0),
      _fid(-1),
      _api_v(0.0),
      _dbs_v(0.0),
      _exErr(0),
      _ifname(std::move(ifname)),
      _isSupported(true),
      _isPopulated(false),
      _isOpen(false),
      _isVerbose(false) {}

exoMesh::~exoMesh() {
  if (_isOpen) ex_close(_fid);
}

void wrnErrMsg(int errCode, const std::string &msg) {
  if (errCode < 0) {
    std::cerr << "Error: " << msg << " with error code " << errCode
              << std::endl;
    throw;
  } else if (errCode > 0) {
    std::cerr << "Warning: " << msg << " with warning code " << errCode
              << std::endl;
  }
}

void exoMesh::write() {
  // preparing database
  // regardless we update it
  exoPopulate(true);

  // writing to file
  int comp_ws = 8;
  int io_ws = 8;
  _fid = ex_create(_ifname.c_str(), EX_CLOBBER, &comp_ws, &io_ws);
  _isOpen = true;

  // initializing exodus database
  _exErr =
      ex_put_init(_fid, "NEMoSys ExodusII database", _numDim, _numNdes,
                  _numElms, _elmBlks.size(), _ndeSets.size(), _sdeSets.size());
  wrnErrMsg(_exErr, "Problem initializing EXODUS II database");

  // writing node coordinates
  _exErr = ex_put_coord(_fid, _xCrds.data(), _yCrds.data(), _zCrds.data());
  wrnErrMsg(_exErr, "Problem writing node coordinates.");

  // writing element blocks
  for (const auto &ieb : _elmBlks) {
    _exErr = ex_put_elem_block(_fid, ieb.id, elmTypeStr(ieb.eTpe).c_str(),
                               ieb.nElm, ieb.ndePerElm, 0);
    wrnErrMsg(_exErr, "Problem writing element block parameters.");
    /*
    int elem_blk_id, num_elem_this_blk, num_nodes_per_elem, num_attr;
    std::string elem_type;
    elem_type.resize(MAX_STR_LENGTH, '\0');
    elem_blk_id = ieb + 1;
    ex_get_elem_block(_fid, elem_blk_id, &elem_type[0], &num_elem_this_blk,
                      &num_nodes_per_elem, &num_attr);
    std::cout << "elem_blk_id : " << elem_blk_id
              << "\nelem_type : " << std::string(elem_type)
              << "\nnum_elem_this_blk : " << num_elem_this_blk
              << "\nnum_nodes_per_elem : " << num_nodes_per_elem
              << "\nnum_attr : " << num_attr
              << std::endl;
    */
    _exErr = ex_put_elem_conn(_fid, ieb.id, ieb.conn.data());
    wrnErrMsg(_exErr, "Problem writing element block connectivities.");
    // write element attribute
    /*
    std::vector<double> elmAttrib;
    elmAttrib.resize(_elmBlock[ieb].nElm, 1);
    _exErr = ex_put_elem_attr(_fid, 1, elmAttrib.data());
    wrnErrMsg(_exErr, "Problem during writing elementBlock attributes.");
    */
  }

  // writing node sets
  for (const auto &ins : _ndeSets) {
    _exErr = ex_put_node_set_param(_fid, ins.id, ins.nNde, 0);
    wrnErrMsg(_exErr, "Problem writing node set parameters.");
    _exErr = ex_put_node_set(_fid, ins.id, ins.ndeIds.data());
    wrnErrMsg(_exErr, "Problem writing node set node ids.");
  }

  // writing side sets
  for (const auto &iss : _sdeSets) {
    _exErr = ex_put_side_set_param(_fid, iss.id, iss.nSde, 0);
    wrnErrMsg(_exErr, "Problem writing side set parameters.");
    _exErr =
        ex_put_side_set(_fid, iss.id, iss.elmIds.data(), iss.sdeIds.data());
    wrnErrMsg(_exErr, "Problem writing side set element and side ids.");
  }

  // writing names
  std::vector<char *> names;
  for (auto &&_elmBlkName : _elmBlkNames) names.push_back(&_elmBlkName[0]);
  if (!names.empty()) {
    _exErr = ex_put_names(_fid, EX_ELEM_BLOCK, names.data());
    wrnErrMsg(_exErr, "Problem writing element block names.");
  }

  names.clear();
  for (auto &&_ndeSetName : _ndeSetNames) names.push_back(&_ndeSetName[0]);
  if (!names.empty()) {
    _exErr = ex_put_names(_fid, EX_NODE_SET, names.data());
    wrnErrMsg(_exErr, "Problem writing node set names.");
  }

  names.clear();
  for (auto &&_sdeSetName : _sdeSetNames) names.push_back(&_sdeSetName[0]);
  if (!names.empty()) {
    _exErr = ex_put_names(_fid, EX_SIDE_SET, names.data());
    wrnErrMsg(_exErr, "Problem writing side set names.");
  }

  // closing the file
  _exErr = ex_close(_fid);
  wrnErrMsg(_exErr, "Problem closing the EXODUS II database.");
  _isOpen = false;
}

void exoMesh::exoPopulate(bool updElmLst) {
  _numElms = 0;
  _elmBlkNames.clear();
  if (updElmLst) {
    for (auto &&ieb : _elmBlks) ieb.elmIds.clear();
    glbConn.clear();
  }
  _ndeSetNames.clear();
  _sdeSetNames.clear();

  // removing empty element blocks
  std::vector<elmBlkType> nebs;
  for (auto &&ieb : _elmBlks)
    if (ieb.nElm > 0) nebs.emplace_back(std::move(ieb));
  _elmBlks = nebs;

  // removing empty node sets
  std::vector<ndeSetType> nnss;
  for (auto &&ins : _ndeSets)
    if (ins.nNde > 0) nnss.emplace_back(std::move(ins));
  _ndeSets = nnss;

  // removing empty side sets
  std::vector<sdeSetType> nsss;
  for (auto &&iss : _sdeSets)
    if (iss.nSde > 0) nsss.emplace_back(std::move(iss));
  _sdeSets = nsss;

  // updating element blocks
  int blkId = 0;  // Reindex the blocks.
  for (auto &&ieb : _elmBlks) {
    if (_reindexBlks) ieb.id = ++blkId;
    _elmBlkNames.emplace_back(ieb.name);

    // offsetting element connectivities
    if (ieb.ndeIdOffset != 0) {
      for (auto &&conn : ieb.conn) conn += ieb.ndeIdOffset;
      ieb.ndeIdOffset = 0;
    }

    // updating global connectivity
    if (updElmLst) {
      ieb.elmIds.clear();
      int nn = ieb.ndePerElm;
      for (int elmIdx = 0; elmIdx < ieb.nElm; ++elmIdx) {
        std::vector<int> elmConn;
        for (int idx = elmIdx * nn; idx < (elmIdx + 1) * nn; ++idx)
          elmConn.emplace_back(ieb.conn[idx]);
        ieb.elmIds.emplace_back(_numElms + elmIdx + 1);
        glbConn.emplace_back(elmConn);
      }

      // sanity check
      if (ieb.nElm != ieb.elmIds.size())
        std::cerr << "WARNING: Malformed element block " << ieb.name << ".\n";
    }

    // sanity check
    if (ieb.nElm * ieb.ndePerElm != ieb.conn.size())
      std::cerr << "WARNING: Malformed element block " << ieb.name << ".\n";

    _numElms += ieb.nElm;
  }

  // updating node sets
  int nsetId = 0;  // Reindex the node sets.
  for (auto &&ins : _ndeSets) {
    ins.id = ++nsetId;
    _ndeSetNames.emplace_back(ins.name);

    // offsetting node ids
    if (ins.ndeIdOffset != 0) {
      for (auto &&it2 : ins.ndeIds) it2 += ins.ndeIdOffset;
      ins.ndeIdOffset = 0;
    }

    // sanity check
    if (ins.nNde != ins.ndeIds.size())
      std::cerr << "WARNING: Malformed node set " << ins.name << ".\n";
  }

  // updating side sets
  int ssetId = 0;  // Reindex the side sets.
  for (auto &&iss : _sdeSets) {
    iss.id = ++ssetId;
    _sdeSetNames.emplace_back(iss.name);

    // offsetting element ids
    if (iss.elmIdOffset != 0) {
      for (auto &&elm : iss.elmIds) elm += iss.elmIdOffset;
      iss.elmIdOffset = 0;
    }

    // sanity check
    if (iss.nSde != iss.sdeIds.size() || iss.nSde != iss.elmIds.size())
      std::cerr << "WARNING: Malformed side set " << iss.name << ".\n";
  }

  if (_isVerbose) report();

  _isPopulated = true;
}

void exoMesh::report() const {
  std::cout << " ----- Exodus II Database Report ----- \n";
  std::cout << "Database: " << _ifname << "\n";
  std::cout << "Dimension: " << _numDim << "\n";
  std::cout << "Nodes: " << _numNdes << "\n";
  std::cout << "Elements: " << _numElms << "\n";
  std::cout << "Element blocks: " << getNumberOfElementBlocks() << "\n";
  std::cout << "Node sets: " << getNumberOfNodeSets() << "\n";
  std::cout << "Side sets: " << getNumberOfSideSets() << "\n";
  if (getNumberOfElementBlocks() > 0) {
    std::cout << "  Blk      nElm     eType               Name\n";
    std::cout << "  ---      ----     -----              ------\n";
    for (const auto &ieb : _elmBlks)
      std::cout << std::setw(5) << ieb.id << std::setw(10) << ieb.nElm
                << std::setw(10) << elmTypeStr(ieb.eTpe) << std::setw(20)
                << ieb.name << "\n";
  }
  if (getNumberOfNodeSets() > 0) {
    std::cout << "  Nde      nNde               Name\n";
    std::cout << "  ---      ----              ------\n";
    for (const auto &ins : _ndeSets)
      std::cout << std::setw(5) << ins.id << std::setw(10) << ins.nNde
                << std::setw(20) << ins.name << "\n";
  }
  if (getNumberOfSideSets() > 0) {
    std::cout << "  Sde      nSde               Name\n";
    std::cout << "  ---      ----              ------\n";
    for (const auto &iss : _sdeSets)
      std::cout << std::setw(5) << iss.id << std::setw(10) << iss.nSde
                << std::setw(20) << iss.name << "\n";
  }
  std::cout << std::flush;
}

//***************** Combine Element Blocks ***************//
//********************************************************//
void exoMesh::combineElmBlks(const std::vector<int> &blkIdLst,
                             const std::string &newName) {
  nemAux::Timer T;
  T.start();
  // Check for duplicate entries in blkIdLst
  if (findDuplicates(blkIdLst).size() != 0) {
    std::cerr << "Error: Duplicate block IDs defined." << std::endl;
    throw;
  }
  // Check there is more than one id in list
  if (blkIdLst.size() < 2) {
    std::cerr << "Error: More than one block ID is needed." << std::endl;
    throw;
  }
  // Check that block ID exists
  elementType elmtype;
  std::vector<elementType> eType_vec;
  for (auto &id : blkIdLst) {
    auto it = std::find(_elmBlkIds.begin(), _elmBlkIds.end(), id);
    if (it == _elmBlkIds.end()) {
      std::cerr << "Error: Block ID " << id << " does not exist." << std::endl;
      throw;
    }
    // Check element types for each ID
    elmtype = getElmBlkTypeById(id);
    eType_vec.push_back(elmtype);
  }
  // Check that element type in each block is the same
  for (int i = 0; i < eType_vec.size(); ++i) {
    elementType e;
    if (i == 0)
      e = eType_vec[i];
    else if (eType_vec[i] != e) {
      std::cerr << "Error: Block element types are not the same. Cannot "
                   "combine blocks."
                << std::endl;
      throw;
    }
  }

  // Loop through registered blocks and make vector list
  std::vector<int> blkIds;
  for (auto &blk : _elmBlks) {
    blkIds.push_back(blk.id);
  }

  std::vector<int> arrange_list;
  // Get the first id in combine block list
  int first = blkIdLst[0];
  // Loop through blkIds and insert into new vector
  arrange_list.push_back(blkIds[0]);
  for (int &id : blkIds) {
    if (id != first) {
      // Check if id is already in arrange_list,
      // if so, don't insert
      auto it = std::find(arrange_list.begin(), arrange_list.end(), id);
      if (it == arrange_list.end()) arrange_list.push_back(id);
    } else {
      // Check if each combine list item is in arrange_list
      // if so, remove from arrange_list
      for (const int &comId : blkIdLst) {
        for (int i = 0; i < arrange_list.size(); ++i) {
          if (comId == arrange_list[i]) {
            arrange_list.erase(arrange_list.begin() + i);
          }
        }
      }
      // Insert entire combine list into arrange_list
      for (const int &j : blkIdLst) {
        arrange_list.push_back(j);
      }
    }
  }

  // Populate old to new element ids map
  std::map<int, int> old2NewElmIds;
  int i = 0;
  int numElms = getNumberOfElements();
  for (int &id : arrange_list) {
    for (auto &elmBlk : _elmBlks) {
      if (id == elmBlk.id) {
        for (auto &eid : elmBlk.elmIds) {
          i++;
          old2NewElmIds.insert(std::pair<int, int>(eid, i));
        }
      }
    }
  }

  //**************************************************************//
  // Any part of a side set that is between the blocks being
  // combined should be deleted. Therefore, delete these parts
  // before combining the blocks by:
  //
  // Find all node IDs common between combining blocks (blkIdLst).
  // Build combining blocks point kdTree.
  // Go block by block and insert points into tree.
  // Find elements that have these common node IDs.
  // Find those elements in the side sets and compare common nodes
  // to the face nodes.
  //**************************************************************//
  std::cout << "\nModifying sidesets..." << std::flush;
  // Get the number nodes in block to be combined
  int numNodes = 0;
  for (auto &id : blkIdLst) {
    int idx = getElmBlkIndex(id);
    numNodes += _elmBlks[idx].ndeCoords.size();
  }

  // Build kdtree from node coords in each block
  // This will create duplicate nodes in the kdtree
  // which is desired to find common nodes b/w blocks
  int nDim = 3;
  std::map<int, int> tree_ndeId_map;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(numNodes, nDim);
  int ipts = 0;
  for (auto &id : blkIdLst) {
    int idx = getElmBlkIndex(id);
    auto it = _elmBlks[idx].ndeCoords.begin();
    while (it != _elmBlks[idx].ndeCoords.end()) {
      pntCrd[ipts][0] = it->second[0];
      pntCrd[ipts][1] = it->second[1];
      pntCrd[ipts][2] = it->second[2];
      tree_ndeId_map[ipts] = it->first;
      ipts++;
      it++;
    }
  }
  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, numNodes, nDim);

  // Finding duplicate node ids
  std::set<int> commonNodes;
  int nNib = 2;
  ANNpoint qryPnt;
  ANNidxArray nnIdx = new ANNidx[nNib];
  ANNdistArray dists = new ANNdist[nNib];
  qryPnt = annAllocPt(nDim);

  for (int iNde = 0; iNde < _numNdes; iNde++) {
    qryPnt[0] = _xCrds[iNde];
    qryPnt[1] = _yCrds[iNde];
    qryPnt[2] = _zCrds[iNde];
    kdTree->annkSearch(qryPnt, nNib, nnIdx, dists);

    if (dists[1] <= 1e-12) {
      int id1 = tree_ndeId_map[nnIdx[0]];
      int id2 = tree_ndeId_map[nnIdx[1]];
      commonNodes.insert(id1);
      commonNodes.insert(id2);
    }
  }
  delete kdTree;

  // Find all elements that have commonNodes
  std::set<int> commonElms;
  int idx, num_nodes_per_elem, ndeID, c;
  numElms = 0;
  for (auto &id : blkIdLst) {
    idx = getElmBlkIndex(id);
    numElms = _elmBlks[idx].nElm;
    num_nodes_per_elem = _elmBlks[idx].ndePerElm;
    c = 0;
    for (int el = 0; el < numElms; ++el) {
      for (int n = 0; n < num_nodes_per_elem; ++n) {
        c = el * num_nodes_per_elem + n;
        ndeID = _elmBlks[idx].conn[c];
        auto pos = commonNodes.find(ndeID);
        if (pos != commonNodes.end())
          commonElms.insert(_elmBlks[idx].elmIds[el]);
      }
    }
  }

  // Face to local element connectivity maps
  std::map<int, std::vector<int>> tri_face_conn_map = {
      {1, {1, 2}}, {2, {2, 3}}, {3, {3, 1}}};
  std::map<int, std::vector<int>> quad_face_conn_map = {
      {1, {1, 2}}, {2, {2, 3}}, {3, {3, 4}}, {4, {4, 1}}};
  std::map<int, std::vector<int>> tet_face_conn_map = {
      {1, {1, 2, 4}}, {2, {2, 3, 4}}, {3, {1, 4, 3}}, {4, {1, 3, 2}}};
  std::map<int, std::vector<int>> hex_face_conn_map = {
      {1, {1, 2, 6, 5}}, {2, {2, 3, 7, 6}}, {3, {3, 4, 8, 7}},
      {4, {1, 5, 8, 4}}, {5, {1, 4, 3, 2}}, {6, {5, 6, 7, 8}}};
  std::map<int, std::vector<int>> wedge_face_conn_map = {{1, {1, 2, 5, 4}},
                                                         {2, {2, 3, 6, 5}},
                                                         {3, {3, 1, 4, 6}},
                                                         {4, {1, 2, 3}},
                                                         {5, {4, 5, 6}}};

  int numCommon = 0;
  // Find commonElms in sides sets
  for (auto &sdeSet : _sdeSets) {
    // Iterators for face id vector
    auto itf = sdeSet.sdeIds.begin();
    int loc = 0;
    int elmID = 0;
    for (auto ite = sdeSet.elmIds.begin(); ite != sdeSet.elmIds.end();) {
      loc = ite - sdeSet.elmIds.begin();
      elmID = sdeSet.elmIds[loc];
      auto pos = commonElms.find(elmID);
      if (pos != commonElms.end()) {
        int faceID = sdeSet.sdeIds[loc];

        // Use the elmID to get the connectivity
        std::vector<int> elmConn = glbConn[elmID - 1];
        std::vector<int> faceConn;
        // With the faceID get the local connectivity (map)
        if (elmtype == TETRA)
          faceConn = tet_face_conn_map[faceID];
        else if (elmtype == HEX)
          faceConn = hex_face_conn_map[faceID];
        else if (elmtype == WEDGE)
          faceConn = wedge_face_conn_map[faceID];
        else if (elmtype == TRIANGLE)
          faceConn = tri_face_conn_map[faceID];
        else if (elmtype == QUAD)
          faceConn = quad_face_conn_map[faceID];
        else {
          std::cerr << "Error: Mesh element type not supported." << std::endl;
          throw;
        }

        // Compare the local(face) nodes with commonNodes
        bool allCommon = false;
        for (int i = 0; i < faceConn.size(); ++i) {
          // This is the node ID on the face
          int ndeID = elmConn[faceConn[i] - 1];
          auto pos = commonNodes.find(ndeID);
          if (pos == commonNodes.end()) {
            // If a commonNode is not found, break
            allCommon = false;
            break;
          } else
            allCommon = true;
        }
        if (allCommon) {
          ite = sdeSet.elmIds.erase(ite);
          itf = sdeSet.sdeIds.erase(itf);
          sdeSet.nSde--;
          numCommon++;
        } else {
          ++ite;
          ++itf;
        }
      } else {
        ++ite;
        ++itf;
      }
    }
  }
  std::cout << "done." << std::endl;

  //**************************************************************//
  // Combine the blocks
  //**************************************************************//
  std::cout << "Combining blocks..." << std::flush;
  int firstBlk = getElmBlkIndex(blkIdLst[0]);
  // Combine the blocks
  for (int i = 0; i < blkIdLst.size(); ++i) {
    if (i > 0) {
      int idx = getElmBlkIndex(blkIdLst[i]);
      int ei = -1;
      for (const auto &eid : _elmBlks[idx].elmIds) {
        ei++;
        _elmBlks[firstBlk].elmIds.emplace_back(eid);
        _elmBlks[firstBlk].nElm++;

        for (int ni = ei * _elmBlks[idx].ndePerElm;
             ni < (ei + 1) * _elmBlks[idx].ndePerElm; ni++) {
          _elmBlks[firstBlk].conn.emplace_back(_elmBlks[idx].conn[ni]);
        }
      }
      // Set new block name, if given
      if (newName != "") _elmBlks[firstBlk].name = newName;
    }
  }
  std::cout << "done." << std::endl;
  std::cout << "Updating sideset element IDs..." << std::flush;
  // Update sides sets
  updateSidesets(old2NewElmIds);
  std::cout << "done." << std::endl;

  // Remove blocks, don't reindex blocks
  _reindexBlks = false;
  for (int i = 0; i < blkIdLst.size(); ++i) {
    if (i > 0) removeElmBlkById(blkIdLst[i]);
  }
  T.stop();

  std::cout << "Combined blocks and updated sidesets in "
            << T.elapsed() / 1000.0 << " seconds.\n"
            << std::endl;
}

//******************** Update Sidesets *******************//
//********************************************************//
void exoMesh::updateSidesets(const std::map<int, int> &old2NewElmIds) {
  // Loop through side sets, update elem ids with map
  for (auto &sSet : _sdeSets) {
    for (auto &elmId : sSet.elmIds) {
      auto it = old2NewElmIds.find(elmId);
      if (it != old2NewElmIds.end())
        elmId = it->second;
      else
        std::cerr << "Error: element not found in sideset" << std::endl;
    }
  }
}

void exoMesh::removeByElmIdLst(int blkIdx, const std::vector<int> &idLst) {
  // TODO: Does NOT update side sets properly!
  // range check
  if (blkIdx > getNumberOfElementBlocks() || blkIdx < 0)
    wrnErrMsg(-1, "Block Id is out of the range.");

  // the idLst should not eliminate the old block
  // at least one element should remain
  // if (idLst.size() >= _elmBlock[blkIdx].elmIds.size())
  //{
  //    std::cout << idLst.size() << " " << _elmBlock[blkIdx].elmIds.size() <<
  //    std::endl; wrnErrMsg(-1, "Can not remove the entire block!.");
  //}

  // now updating the old block
  _isPopulated = false;

  elmBlkType oeb = _elmBlks[blkIdx];
  elmBlkType neb;
  neb.id = oeb.id;
  neb.name = oeb.name;
  neb.eTpe = oeb.eTpe;
  neb.ndePerElm = oeb.ndePerElm;
  neb.ndeIdOffset = 0;
  neb.nElm = 0;

  // elmIds and connectivities
  int oei = -1;
  for (const auto &eid : oeb.elmIds) {
    oei++;

    bool stays = true;
    for (const auto &id : idLst)
      if (id == eid) {
        stays = false;
        break;
      }
    if (!stays) continue;

    neb.nElm++;
    neb.elmIds.emplace_back(eid);
    for (int ni = oei * oeb.ndePerElm; ni < (oei + 1) * oeb.ndePerElm; ni++)
      neb.conn.emplace_back(oeb.conn[ni]);
  }
  std::cout << "Removed " << oeb.elmIds.size() - neb.elmIds.size()
            << " elements from original" << std::endl;
  _elmBlks[blkIdx] = neb;
}

void exoMesh::addElmBlkByElmIdLst(const std::string &name,
                                  std::vector<int> &lst) {
  if (lst.empty()) {
    std::cerr << "WARNING: Attempting to add element block " << name
              << " with no elements." << std::endl;
    return;
  }

  // preparing database
  if (!_isPopulated) exoPopulate(false);

  _isPopulated = false;
  // create element block
  // Assumption here: All elements in the provided list are in the same element
  // block
  elmBlkType neb;
  neb.id = _elmBlks.size() + 1;
  neb.name = name;
  neb.ndeIdOffset = 0;

  // int blkIdx = findElmBlkIdxByElmId(lst[0]); // faster but less accurate
  int blkIdx = findElmBlkIdxByElmIdLst(lst);  // slower, more accurate
  if (blkIdx == -1) wrnErrMsg(-1, "Elements ids are not registered.");
  // std::cout << "Block Indx = " << blkIdx << std::endl;

  // now adjust the list and remove elements that are not owned by this block
  std::vector<int> owndLst;
  bool allIn = false;
  owndLst = lstElmInBlk(blkIdx, lst, allIn);
  if (!allIn) {
    std::cout << owndLst.size() << " Elements are in the block.\n";
    std::cout << "Some of the elements in the list are outside the block.\n";
    lst = owndLst;
  }

  neb.elmIds = lst;
  neb.nElm = lst.size();
  neb.eTpe = getElmBlkType(blkIdx);
  neb.ndePerElm = _elmBlks[blkIdx].ndePerElm;
  for (const auto &eid : lst)
    for (int idx = 0; idx < neb.ndePerElm; ++idx)
      neb.conn.emplace_back(glbConn[eid][idx]);

  // updating the old element block and adding a new one
  removeByElmIdLst(blkIdx, lst);
  addElmBlk(neb);

  exoPopulate(false);
}

void exoMesh::addNdeSetByNdeIdLst(const std::string &name,
                                  const std::vector<int> &idLst) {
  // preparing database
  if (!_isPopulated) exoPopulate(false);

  _isPopulated = false;
  // create and add new node set
  ndeSetType nns;
  nns.id = _ndeSets.size() + 1;  // 1-indexed
  nns.nNde = idLst.size();
  nns.name = name;
  nns.ndeIdOffset = 0;
  nns.ndeIds = idLst;
  addNdeSet(nns);

  // update exodus information
  exoPopulate(false);
}

void exoMesh::snapNdeCrdsZero(double tol) {
  int nPnt = 0;

  // loop through nodes
  bool chk;
  for (auto i = 0; i < _numNdes; ++i) {
    chk = false;
    if (std::abs(_xCrds[i]) <= tol) {
      _xCrds[i] = 0.0;
      chk = true;
    }
    if (std::abs(_yCrds[i]) <= tol) {
      _yCrds[i] = 0.0;
      chk = true;
    }
    if (std::abs(_zCrds[i]) <= tol) {
      _zCrds[i] = 0.0;
      chk = true;
    }
    if (chk) nPnt++;
  }

  std::cout << "Number of points snapped : " << nPnt << std::endl;
}

void exoMesh::removeElmBlkByName(const std::string &blkName) {
  if (blkName.empty()) return;

  std::vector<elmBlkType> nebs;
  for (const auto &ieb : _elmBlks)
    if (ieb.name != blkName)
      nebs.emplace_back(ieb);
    else
      _isPopulated = false;
  _elmBlks = nebs;

  exoPopulate(true);
}

void exoMesh::removeElmBlkById(int id) {
  std::vector<elmBlkType> nebs;
  for (const auto &ieb : _elmBlks)
    if (ieb.id != id)
      nebs.emplace_back(ieb);
    else
      _isPopulated = false;
  _elmBlks = nebs;

  exoPopulate(true);
}

int exoMesh::findElmBlkIdxByElmId(int elmId) const {
  int blkIdx = 0;
  for (const auto &ieb : _elmBlks) {
    for (const auto &eid : ieb.elmIds)
      if (eid == elmId) return blkIdx;

    blkIdx++;
  }
  return -1;
}

int exoMesh::findElmBlkIdxByElmIdLst(const std::vector<int> &elmIds) const {
  if (elmIds.empty()) return -1;

  int blkIdx = -1;
  int intfCnt = 0;

  bool allIn = false;
  for (int ieb = 0; ieb < _elmBlks.size(); ++ieb) {
    std::vector<int> intfLst;
    intfLst = lstElmInBlk(ieb, elmIds, allIn);

    if (allIn) {
      return ieb;
    } else if (intfLst.size() > intfCnt) {
      intfCnt = intfLst.size();
      blkIdx = ieb;
    }
  }

  if (!allIn)
    std::cerr << "WARNING: Elements in list are not all in one element block. "
                 "Only block containing largest amount of elements is "
                 "processed. This might indicate selected elements overlap or "
                 "are different element types (e.g., HEX and TETRA)."
              << std::endl;

  return blkIdx;
}

std::vector<int> exoMesh::lstElmInBlk(int blkIdx,
                                      const std::vector<int> &elmIds,
                                      bool &allIn) const {
  std::vector<int> out;
  allIn = true;

  // compare elmIds with list of the elements within block
  std::set<int> current;
  for (const auto &elmId : _elmBlks[blkIdx].elmIds) current.insert(elmId);
  for (const auto &elmId : elmIds)
    if (!current.insert(elmId).second)
      out.emplace_back(elmId);
    else
      allIn = false;

  return out;
}

void exoMesh::reset() {
  _isPopulated = false;
  _ndeSets.clear();
  _elmBlks.clear();
  _sdeSets.clear();
  _fid = 0;
  _numDim = 2;
  _numNdes = 0;
  _numElms = 0;
  _xCrds.clear();
  _yCrds.clear();
  _zCrds.clear();
  glbConn.clear();
  _ndeSetNames.clear();
  _elmBlkNames.clear();
  _sdeSetNames.clear();
}

void exoMesh::read(const std::string &ifname) {
  if (!ifname.empty()) _ifname = ifname;

  // before reading all internal data base will be reset
  reset();

  int CPU_word_size = sizeof(double);
  int IO_word_size = 0;
  float version;
  _exErr = 0;

  /* open EXODUS II files */
  _fid = ex_open(_ifname.c_str(), EX_READ, &CPU_word_size, &IO_word_size,
                 &version);
  _isOpen = true;
  wrnErrMsg(_fid > 0 ? 0 : -1, "Problem opening file " + _ifname + "\n");

  int numElmBlks;
  int numNdeSets;
  int numSdeSets;

  // parameter inquiry from Exodus file
  int idum;
  float fdum;
  char cdum;

  _exErr = ex_inquire(_fid, EX_INQ_API_VERS, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II API version is " << fdum << std::endl;
  _api_v = fdum;

  _exErr = ex_inquire(_fid, EX_INQ_DB_VERS, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II Database version is " << fdum << std::endl;
  _dbs_v = fdum;

  _exErr = ex_inquire(_fid, EX_INQ_DIM, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  _numDim = idum;
  std::cout << "Number of coordinate dimensions is " << idum << std::endl;

  _exErr = ex_inquire(_fid, EX_INQ_NODES, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of points " << idum << std::endl;
  _numNdes = idum;

  _exErr = ex_inquire(_fid, EX_INQ_ELEM, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of elements " << idum << std::endl;
  _numElms = idum;

  _exErr = ex_inquire(_fid, EX_INQ_ELEM_BLK, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of element blocks " << idum << std::endl;
  numElmBlks = idum;

  _exErr = ex_inquire(_fid, EX_INQ_NODE_SETS, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of node sets " << idum << std::endl;
  numNdeSets = idum;

  _exErr = ex_inquire(_fid, EX_INQ_SIDE_SETS, &idum, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of side sets " << idum << std::endl;
  numSdeSets = idum;

  // nodal coordinates
  _xCrds.resize(_numNdes, 0.0);
  _yCrds.resize(_numNdes, 0.0);
  _zCrds.resize(_numNdes, 0.0);
  _exErr = ex_get_coord(_fid, _xCrds.data(), _yCrds.data(), _zCrds.data());
  wrnErrMsg(_exErr, "Problem reading nodal coordinates.\n");

  // ids
  std::vector<int> elmBlkIds(numElmBlks, 0);
  if (numElmBlks > 0) {
    _exErr = ex_get_elem_blk_ids(_fid, elmBlkIds.data());
    wrnErrMsg(_exErr, "Problem reading element block ids.\n");

    for (int i = 0; i < numElmBlks; ++i) {
      _elmBlkIds.push_back(elmBlkIds[i]);
    }
  }

  std::vector<int> ndeSetIds(numNdeSets, 0);
  if (numNdeSets > 0) {
    _exErr = ex_get_node_set_ids(_fid, ndeSetIds.data());
    wrnErrMsg(_exErr, "Problem reading node set ids.\n");
  }

  std::vector<int> sdeSetIds(numSdeSets, 0);
  if (numSdeSets > 0) {
    _exErr = ex_get_side_set_ids(_fid, sdeSetIds.data());
    wrnErrMsg(_exErr, "Problem reading side set ids.\n");
  }

  // element block names
  std::string blk_name;
  for (int iEB = 0; iEB < numElmBlks; ++iEB) {
    blk_name.clear();
    blk_name.resize(MAX_STR_LENGTH);
    int blk_id = elmBlkIds[iEB];
    _exErr = ex_get_name(_fid, EX_ELEM_BLOCK, blk_id, &blk_name[0]);
    wrnErrMsg(_exErr, "Problem reading element block names.\n");
    blk_name.erase(std::remove(blk_name.begin(), blk_name.end(), '\0'),
                   blk_name.end());
    // std::cout << "debugging" << std::endl;
    // std::cout << blk_id << " " << blk_name << std::endl;
    _elmBlkNames.push_back(blk_name);
  }

  // node set names
  for (int iNS = 0; iNS < numNdeSets; ++iNS) {
    blk_name.clear();
    blk_name.resize(MAX_STR_LENGTH);
    int blk_id = ndeSetIds[iNS];
    _exErr = ex_get_name(_fid, EX_NODE_SET, blk_id, &blk_name[0]);
    wrnErrMsg(_exErr, "Problem reading node set names.\n");
    blk_name.erase(std::remove(blk_name.begin(), blk_name.end(), '\0'),
                   blk_name.end());
    // std::cout << blk_id << " " << blk_name << std::endl;
    _ndeSetNames.push_back(blk_name);
  }

  // side set names
  for (int iSS = 0; iSS < numSdeSets; ++iSS) {
    blk_name.clear();
    blk_name.resize(MAX_STR_LENGTH);
    int blk_id = sdeSetIds[iSS];
    _exErr = ex_get_name(_fid, EX_SIDE_SET, blk_id, &blk_name[0]);
    wrnErrMsg(_exErr, "Problem reading side set names.\n");
    blk_name.erase(std::remove(blk_name.begin(), blk_name.end(), '\0'),
                   blk_name.end());
    // std::cout << blk_id << " " << blk_name << std::endl;
    _sdeSetNames.push_back(blk_name);
  }

  // element blocks
  std::string elem_type;
  int elmId = 1;
  for (int iEB = 0; iEB < numElmBlks; ++iEB) {
    elem_type.clear();
    elem_type.resize(MAX_STR_LENGTH);
    elmBlkType eb;
    eb.id = elmBlkIds[iEB];
    int num_el_in_blk, num_nod_per_el, num_attr;
    _exErr = ex_get_elem_block(_fid, elmBlkIds[iEB], &elem_type[0],
                               &num_el_in_blk, &num_nod_per_el, &num_attr);
    wrnErrMsg(_exErr, "Problem reading element block parameters.\n");
    elem_type.erase(std::remove(elem_type.begin(), elem_type.end(), '\0'),
                    elem_type.end());
    eb.nElm = num_el_in_blk;
    eb.name = _elmBlkNames[iEB];
    eb.ndePerElm = num_nod_per_el;
    eb.eTpe = elmTypeNum(elem_type);

    eb.ndeIdOffset = 0;

    // read element connectivity
    eb.conn.clear();
    eb.conn.resize(num_el_in_blk * num_nod_per_el, 0);
    _exErr = ex_get_elem_conn(_fid, elmBlkIds[iEB], eb.conn.data());
    wrnErrMsg(_exErr, "Problem reading element block connectivities.\n");

    double x, y, z;
    // get node coords for each block
    for (int &ndeId : eb.conn) {
      x = _xCrds[ndeId];
      y = _yCrds[ndeId];
      z = _zCrds[ndeId];
      eb.ndeCoords.insert(
          std::pair<int, std::vector<double>>(ndeId, {x, y, z}));
    }

    // fill element ids
    eb.elmIds.clear();
    eb.elmIds.resize(num_el_in_blk);
    std::iota(eb.elmIds.begin(), eb.elmIds.end(), elmId);
    elmId += num_el_in_blk;

    // global connectivity
    int nn = eb.ndePerElm;
    for (int elmIdx = 0; elmIdx < eb.nElm; ++elmIdx) {
      std::vector<int> elmConn;
      for (int idx = elmIdx * nn; idx < (elmIdx + 1) * nn; ++idx)
        elmConn.emplace_back(eb.conn[idx]);
      glbConn.emplace_back(elmConn);
    }

    _elmBlks.push_back(eb);
  }

  // node sets
  for (int iNS = 0; iNS < numNdeSets; ++iNS) {
    ndeSetType ns;
    ns.id = ndeSetIds[iNS];
    _exErr = ex_get_node_set_param(_fid, ns.id, &ns.nNde, &idum);
    wrnErrMsg(_exErr, "Problem reading node set parameters.\n");
    ns.name = _ndeSetNames[iNS];

    ns.ndeIdOffset = 0;
    ns.ndeIds.resize(ns.nNde, 0);
    _exErr = ex_get_node_set(_fid, ns.id, ns.ndeIds.data());
    wrnErrMsg(_exErr, "Problem reading node set node ids.\n");
    _ndeSets.push_back(ns);
  }

  // side sets
  for (int iSS = 0; iSS < numSdeSets; ++iSS) {
    sdeSetType ss;
    ss.id = sdeSetIds[iSS];
    _exErr = ex_get_side_set_param(_fid, ss.id, &ss.nSde, &idum);
    wrnErrMsg(_exErr, "Problem reading side set parameters.\n");
    ss.name = _sdeSetNames[iSS];

    ss.elmIdOffset = 0;
    ss.sdeIds.resize(ss.nSde, 0);
    ss.elmIds.resize(ss.nSde, 0);
    _exErr = ex_get_side_set(_fid, ss.id, ss.elmIds.data(), ss.sdeIds.data());
    wrnErrMsg(_exErr, "Problem reading side set element and side ids.\n");
    _sdeSets.push_back(ss);
  }

  report();

  // since all data structures are manually populated from the file
  _isPopulated = true;
}

void exoMesh::mergeNodes(double tol) {
  // build all point kdTree
  int nDim = 3;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(_numNdes, nDim);
  for (int iPnt = 0; iPnt < _numNdes; iPnt++) {
    pntCrd[iPnt][0] = _xCrds[iPnt];
    pntCrd[iPnt][1] = _yCrds[iPnt];
    pntCrd[iPnt][2] = _zCrds[iPnt];
  }
  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, _numNdes, nDim);

  // finding duplicate node ids
  std::map<int, int> dupNdeMap;  // exoMesh is one-indexed
  int nNib = 2;
  ANNpoint qryPnt;
  ANNidxArray nnIdx = new ANNidx[nNib];
  ANNdistArray dists = new ANNdist[nNib];
  qryPnt = annAllocPt(nDim);
  for (int iNde = 0; iNde < _numNdes; iNde++) {
    qryPnt[0] = _xCrds[iNde];
    qryPnt[1] = _yCrds[iNde];
    qryPnt[2] = _zCrds[iNde];
    kdTree->annkSearch(qryPnt, nNib, nnIdx, dists);
    if (dists[1] <= tol) dupNdeMap[nnIdx[0] + 1] = nnIdx[1] + 1;
  }
  std::cout << "Found " << dupNdeMap.size() << " duplicate nodes.\n";

  // removing duplicated nodal coordinates
  std::vector<double> xn, yn, zn;
  std::map<int, int> old2NewNde;
  for (int iNde = 0; iNde < _numNdes; iNde++) {
    if (dupNdeMap.find(iNde + 1) == dupNdeMap.end()) {
      xn.push_back(_xCrds[iNde]);
      yn.push_back(_yCrds[iNde]);
      zn.push_back(_zCrds[iNde]);
    }
  }

  delete kdTree;

  // mapping
  for (int iPnt = 0; iPnt < xn.size(); iPnt++) {
    pntCrd[iPnt][0] = xn[iPnt];
    pntCrd[iPnt][1] = yn[iPnt];
    pntCrd[iPnt][2] = zn[iPnt];
  }
  kdTree = new ANNkd_tree(pntCrd, xn.size(), nDim);

  // creating map from old ids to new ids
  for (int iNde = 0; iNde < _numNdes; iNde++) {
    qryPnt[0] = _xCrds[iNde];
    qryPnt[1] = _yCrds[iNde];
    qryPnt[2] = _zCrds[iNde];
    kdTree->annkSearch(qryPnt, 1, nnIdx, dists);
    if (dists[0] <= tol) {
      old2NewNde[iNde + 1] = nnIdx[0] + 1;
    } else {
      std::cerr << "ERROR: Found a node in database not properly copied.\n";
      exit(1);
    }
  }

  // sanity check
  int max = 0;
  int min = std::numeric_limits<int>::max();
  for (const auto &im : old2NewNde) {
    max = std::max(max, im.second);
    min = std::min(min, im.second);
  }
  std::cout << "Max new node id = " << max << "\n";
  std::cout << "Min new node id = " << min << "\n";

  // update node coordinates
  _xCrds.clear();
  _yCrds.clear();
  _zCrds.clear();
  _xCrds.assign(xn.begin(), xn.end());
  _yCrds.assign(yn.begin(), yn.end());
  _zCrds.assign(zn.begin(), zn.end());
  std::cout << "Number of nodes changed from " << _numNdes << " to "
            << xn.size() << "\n";
  _numNdes = xn.size();

  // update element block connectivities
  max = 0;
  min = std::numeric_limits<int>::max();
  for (auto &&ieb : _elmBlks) {
    for (auto &&icn : ieb.conn) {
      int nid = icn;
      icn = old2NewNde[nid];
      max = std::max(max, icn);
      min = std::min(min, icn);
    }
  }
  std::cout << "Max conn = " << max << "\n";
  std::cout << "Min conn = " << min << "\n";

  // update node set ids
  for (auto &&ins : _ndeSets) {
    // updating node indices, removing duplicates
    std::set<int> newNSNdeIds;  // set must have unique entries

    for (const auto &nid : ins.ndeIds) newNSNdeIds.insert(old2NewNde[nid]);

    ins.ndeIds.clear();
    ins.ndeIds.assign(newNSNdeIds.begin(), newNSNdeIds.end());
    std::cout << "Original node set with " << ins.nNde
              << " nodes after cleaning becomes " << ins.ndeIds.size() << "\n";

    // updating the nodeset data structure
    ins.nNde = ins.ndeIds.size();
  }

  // side sets do not need to be updated since they contain only element ids

  // deleting KDTree
  delete kdTree;
  delete[] nnIdx;
  delete[] dists;

  // Need to update the element list.
  exoPopulate(true);
}

void exoMesh::scaleNodes(double sc) {
  using namespace nemAux;
  _xCrds = sc * _xCrds;
  _yCrds = sc * _yCrds;
  _zCrds = sc * _zCrds;
}

void exoMesh::stitch(const exoMesh &otherMesh) {
  // Append nodes.
  _xCrds.insert(_xCrds.end(), otherMesh._xCrds.begin(), otherMesh._xCrds.end());
  _yCrds.insert(_yCrds.end(), otherMesh._yCrds.begin(), otherMesh._yCrds.end());
  _zCrds.insert(_zCrds.end(), otherMesh._zCrds.begin(), otherMesh._zCrds.end());

  // Append all element blocks with extra offset.
  for (const auto &elmBlk : otherMesh._elmBlks) {
    // offset by number of nodes in current mesh.
    auto offElmBlk = elmBlk;
    offElmBlk.ndeIdOffset += _numNdes;

    addElmBlk(offElmBlk);
  }

  // Append all node sets.
  for (const auto &ndeSet : otherMesh._ndeSets) {
    auto offNdeSet = ndeSet;
    offNdeSet.ndeIdOffset += _numNdes;

    addNdeSet(offNdeSet);
  }

  // Append all side sets.
  for (const auto &sdeSet : otherMesh._sdeSets) {
    auto offSdeSet = sdeSet;
    offSdeSet.elmIdOffset += _numElms;

    addSdeSet(offSdeSet);
  }

  // Increase node count. This is delayed as the old value is needed for adding
  // offset to blocks and sets.
  _numNdes += otherMesh._numNdes;

  // Full update to the database. Also updates element count _numElms
  exoPopulate(true);
}

elementType exoMesh::getElmBlkTypeById(int id) const {
  auto eb = _elmBlks.begin();
  for (; eb != _elmBlks.end(); eb++)
    if (eb->id == id) return (eb->eTpe);

  // warning if the element block id is not registered
  if (eb == _elmBlks.end())
    std::cerr << "Warning: There is no element block with id " << id
              << std::endl;

  return (NEM::MSH::EXOMesh::elementType::OTHER);
}

elementType exoMesh::getElmBlkType(const std::string &ebName) const {
  auto eb = _elmBlks.begin();
  for (; eb != _elmBlks.end(); eb++)
    if (eb->name == ebName) return (eb->eTpe);

  // warning if the element block name is not registered
  if (eb == _elmBlks.end())
    std::cerr << "Warning: There is no element block with name " << ebName
              << std::endl;

  return (NEM::MSH::EXOMesh::elementType::OTHER);
}

int exoMesh::getNumElmsInBlk(const std::string &ebName) const {
  auto eb = _elmBlks.begin();
  for (; eb != _elmBlks.end(); eb++)
    if (eb->name == ebName) return (eb->nElm);

  // warning if the element block name is not registered
  if (eb == _elmBlks.end())
    std::cerr << "Warning: There is no element block with name " << ebName
              << std::endl;

  return -1;
}

int exoMesh::getNumElmsInBlkById(int id) const {
  auto eb = _elmBlks.begin();
  for (; eb != _elmBlks.end(); eb++)
    if (eb->id == id) return (eb->nElm);

  // warning if the element block name is not registered
  if (eb == _elmBlks.end())
    std::cerr << "Warning: There is no element block with ID " << id
              << std::endl;

  return -1;
}

const std::string exoMesh::getElmBlkNameById(int id) const {
  auto eb = _elmBlks.begin();
  for (; eb != _elmBlks.end(); eb++)
    if (eb->id == id) return (eb->name);

  // warning if the element block name is not registered
  if (eb == _elmBlks.end())
    std::cerr << "Warning: There is no element block with ID " << id
              << std::endl;

  return "";
}

int exoMesh::getElmBlkIndex(int id) const {
  int numElmBlks = getNumberOfElementBlocks();
  for (int i = 0; i < numElmBlks; ++i) {
    if (_elmBlks[i].id == id) return i;
  }
}

int exoMesh::getElmBlkIndex(std::string name) const {
  int numElmBlks = getNumberOfElementBlocks();
  for (int i = 0; i < numElmBlks; ++i) {
    if (_elmBlks[i].name == name) return i;
  }
}

int exoMesh::getElmBlkId(std::string ebName) const {
  auto eb = _elmBlks.begin();
  for (; eb != _elmBlks.end(); eb++)
    if (eb->name == ebName) return (eb->id);

  // warning if the element block name is not registered
  if (eb == _elmBlks.end())
    std::cerr << "Warning: There is no element block with name " << ebName
              << std::endl;

  return -1;
}

std::string exoMesh::getNdeSetNameById(int id) const {
  auto ns = _ndeSets.begin();
  for (; ns != _ndeSets.end(); ns++)
    if (ns->id == id) return (ns->name);

  // warning if the element block name is not registered
  if (ns == _ndeSets.end())
    std::cerr << "Warning: There is no nodeset with ID " << id << std::endl;

  return "";
}

int exoMesh::getNumNdesInNdeSet(const std::string &nsName) const {
  auto ns = _ndeSets.begin();
  for (; ns != _ndeSets.end(); ns++)
    if (ns->name == nsName) return (ns->nNde);

  // warning if the element block name is not registered
  if (ns == _ndeSets.end())
    std::cerr << "Warning: There is no nodeset with name " << nsName
              << std::endl;

  return -1;
}

int exoMesh::getNumNdesInNdeSetById(int id) const {
  auto ns = _ndeSets.begin();
  for (; ns != _ndeSets.end(); ns++)
    if (ns->id == id) return (ns->nNde);

  // warning if the element block name is not registered
  if (ns == _ndeSets.end())
    std::cerr << "Warning: There is no nodeset with ID " << id << std::endl;

  return -1;
}

int exoMesh::getNdeSetId(const std::string &nsName) const {
  auto ns = _ndeSets.begin();
  for (; ns != _ndeSets.end(); ns++)
    if (ns->name == nsName) return (ns->id);

  // warning if the element block name is not registered
  if (ns == _ndeSets.end())
    std::cerr << "Warning: There is no nodeset with name " << nsName
              << std::endl;

  return -1;
}

int exoMesh::getNdeSetIndex(const std::string &nsName) const {
  auto ns = _ndeSets.begin();
  for (; ns != _ndeSets.end(); ns++)
    if (ns->name == nsName) return ns - _ndeSets.begin();

  // warning if the element block name is not registered
  if (ns == _ndeSets.end())
    std::cerr << "Warning: There is no nodeset with name " << nsName
              << std::endl;

  return -1;
}

int exoMesh::getSdeSetId(const std::string &ssName) const {
  auto ss = _sdeSets.begin();
  for (; ss != _sdeSets.end(); ss++)
    if (ss->name == ssName) return (ss->id);

  // warning if the sideset name is not registered
  if (ss == _sdeSets.end())
    std::cerr << "Warning: There is no sideset with name " << ssName
              << std::endl;

  return -1;
}

int exoMesh::getSdeSetIndex(const std::string &ssName) const {
  auto ss = _sdeSets.begin();
  for (; ss != _sdeSets.end(); ss++)
    if (ss->name == ssName) return ss - _sdeSets.begin();

  // warning if the sideset name is not registered
  if (ss == _sdeSets.end())
    std::cerr << "Warning: There is no sideset with name " << ssName
              << std::endl;

  return -1;
}

int exoMesh::getNumSdesInSdeSet(const std::string &ssName) const {
  auto ss = _sdeSets.begin();
  for (; ss != _sdeSets.end(); ss++)
    if (ss->name == ssName) return (ss->nSde);

  // warning if the element block name is not registered
  if (ss == _sdeSets.end())
    std::cerr << "Warning: There is no sideset with name " << ssName
              << std::endl;

  return -1;
}

int exoMesh::getNumSdesInSdeSetById(int id) const {
  auto ss = _sdeSets.begin();
  for (; ss != _sdeSets.end(); ss++)
    if (ss->id == id) return (ss->nSde);

  // warning if the element block name is not registered
  if (ss == _sdeSets.end())
    std::cerr << "Warning: There is no sideset with ID " << id << std::endl;

  return -1;
}

}  // namespace EXOMesh
}  // namespace MSH
}  // namespace NEM
