#ifdef HAVE_EXODUSII
#include "AuxiliaryFunctions.H"
#include "meshBase.H"
#include "exodusII.h"
#include "exoMesh.H"

// third party
#include <ANN/ANN.h>

#include <set>
#include <fstream>
#include <algorithm>
#include <math.h>

using namespace EXOMesh;

VTKCellType EXOMesh::e2vEMap(elementType et)
{
  std::map<elementType, VTKCellType> eMap =
      {
          {elementType::QUAD,     VTK_QUAD},
          {elementType::TRIANGLE, VTK_TRIANGLE},
          {elementType::TETRA,    VTK_TETRA},
          {elementType::WEDGE,    VTK_WEDGE},
          {elementType::HEX,      VTK_HEXAHEDRON},
          {elementType::OTHER,    VTK_EMPTY_CELL}
      };
  return eMap[et];
}

elementType EXOMesh::v2eEMap(VTKCellType vt)
{
  std::map<VTKCellType, elementType> eMap =
      {
          {VTK_QUAD,                 elementType::QUAD},
          {VTK_QUADRATIC_QUAD,       elementType::QUAD},
          {VTK_TRIANGLE,             elementType::TRIANGLE},
          {VTK_QUADRATIC_TRIANGLE,   elementType::TRIANGLE},
          {VTK_HEXAHEDRON,           elementType::HEX},
          {VTK_QUADRATIC_HEXAHEDRON, elementType::HEX},
          {VTK_TETRA,                elementType::TETRA},
          {VTK_QUADRATIC_TETRA,      elementType::TETRA},
          {VTK_WEDGE,                elementType::WEDGE},
          {VTK_EMPTY_CELL,           elementType::OTHER}
      };
  return eMap[vt];
}

surfaceBCTag EXOMesh::bcTagNum(std::string &tag)
{
  nemAux::toLower(tag);
  if (tag == "fixed") return surfaceBCTag::FIXED;
  if (tag == "symmx") return surfaceBCTag::SYMMX;
  if (tag == "symmy") return surfaceBCTag::SYMMY;
  if (tag == "symmz") return surfaceBCTag::SYMMZ;
  std::cerr << "Unknown surface tag " << tag << std::endl;
  throw;
}

std::string EXOMesh::bcTagStr(int tag)
{
  if (tag == FIXED) return "FIXED";
  if (tag == SYMMX) return "SYMMX";
  if (tag == SYMMY) return "SYMMY";
  if (tag == SYMMZ) return "SYMMZ";
  std::cerr << "Unknown surface tag " << tag << std::endl;
  throw;
}

elementType EXOMesh::elmTypeNum(std::string tag)
{
  nemAux::toLower(tag);
  if (tag == "quadrilateral") return elementType::QUAD;
  if (tag == "quad") return elementType::QUAD;
  if (tag == "triangle") return elementType::TRIANGLE;
  if (tag == "tri") return elementType::TRIANGLE;
  if (tag == "hexahedron") return elementType::HEX;
  if (tag == "hexahedral") return elementType::HEX;
  if (tag == "brick") return elementType::HEX;
  if (tag == "hex") return elementType::HEX;
  if (tag == "tetrahedron") return elementType::TETRA;
  if (tag == "tetrahedral") return elementType::TETRA;
  if (tag == "tetra") return elementType::TETRA;
  if (tag == "wedge") return elementType::WEDGE;
  if (tag == "prismatic") return elementType::WEDGE;
  if (tag == "prism") return elementType::WEDGE;
  std::cout
      << "Warning : Element type "
      << tag
      << " may be not supported\n";
  return elementType::OTHER;
}

std::string EXOMesh::elmTypeStr(elementType tag)
{
  if (tag == elementType::QUAD) return "QUAD";
  if (tag == elementType::TRIANGLE) return "TRIANGLE";
  if (tag == elementType::HEX) return "HEX";
  if (tag == elementType::TETRA) return "TETRA";
  if (tag == elementType::WEDGE) return "WEDGE";
  return "OTHER";
}

int EXOMesh::elmNumNde(elementType tag, int order)
{
  if (tag == elementType::QUAD && order == 1) return 4;
  if (tag == elementType::QUAD && order == 2) return 9;
  if (tag == elementType::TRIANGLE && order == 1) return 3;
  if (tag == elementType::TRIANGLE && order == 2) return 6;
  if (tag == elementType::HEX && order == 1) return 8;
  if (tag == elementType::HEX && order == 2) return 21;
  if (tag == elementType::TETRA && order == 1) return 4;
  if (tag == elementType::TETRA && order == 2) return 10;
  if (tag == elementType::WEDGE && order == 1) return 6;
  if (tag == elementType::WEDGE && order == 2) return 16;
  std::cerr << "Unknown element/order combination\n";
  throw;
}

int EXOMesh::elmNumSrf(elementType tag)
{
  if (tag == elementType::QUAD) return 4;
  if (tag == elementType::TRIANGLE) return 3;
  if (tag == elementType::HEX) return 8;
  if (tag == elementType::TETRA) return 4;
  if (tag == elementType::WEDGE) return 5;
  std::cerr << "Unknown element type\n";
  throw;
}

//////////////////////////////////
// exoMesh class 
//////////////////////////////////

exoMesh::exoMesh(std::string ifname) :
    _ifname(ifname), _isSupported(true), _isPopulated(false),
    _isOpen(true), _isVerbose(false), _avoidPopulation(false),
    _exErr(0)
{}

exoMesh::~exoMesh()
{
  if (_isOpen)
    ex_close(fid);
}


void EXOMesh::wrnErrMsg(int errCode, std::string msg)
{
  if (errCode < 0) {
    std::cerr << "Error: " << msg
              << " with error code "
              << errCode << std::endl;
    throw;
  } else if (errCode > 0) {
    std::cerr << "Warning: code number " << errCode << std::endl;
  }
}

void exoMesh::write()
{
  // preparing database
  // regardless we update it
  exoPopulate(true);
  std::cout << "nNodes = " << numNdes << std::endl;

  // writing to file
  int comp_ws = 8;
  int io_ws = 8;
  fid = ex_create(_ifname.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

  // initializing exodus database
  _exErr = ex_put_init(fid, "NEMoSys ExodusII database", 3, numNdes, numElms,
                       _elmBlock.size(), _ndeSet.size(), _sideSet.size());
  wrnErrMsg(_exErr, "Problem during initialization of exodusII database");

  // writing node coordinates
  _exErr = ex_put_coord(fid, &xCrds[0], &yCrds[0], &zCrds[0]);
  wrnErrMsg(_exErr, "Problem during writing node coordinates.");

  // writing node sets
  for (int ins = 0; ins < _ndeSet.size(); ins++) {
    _exErr = ex_put_node_set_param(fid, _ndeSet[ins].id, _ndeSet[ins].nNde, 0);
    wrnErrMsg(_exErr, "Problem during writing nodeSet parameters.");
    _exErr = ex_put_node_set(fid, _ndeSet[ins].id, &(_ndeSet[ins].ndeIds[0]));
    wrnErrMsg(_exErr, "Problem during writing nodeSet node Ids.");
  }

  // writing element blocks
  for (int ieb = 0; ieb < _elmBlock.size(); ieb++) {
    std::string eTpe = elmTypeStr(_elmBlock[ieb].eTpe);
    _exErr = ex_put_elem_block(fid, _elmBlock[ieb].id, eTpe.c_str(),
                               _elmBlock[ieb].nElm, _elmBlock[ieb].ndePerElm,
                               1);
    wrnErrMsg(_exErr, "Problem during writing elementBlock parameters.");
    int elem_blk_id, num_elem_this_blk, num_nodes_per_elem, num_attr;
    std::string elem_type;
    elem_type.resize(MAX_STR_LENGTH, '\0');
    elem_blk_id = ieb + 1;
    ex_get_elem_block(fid, elem_blk_id, &elem_type[0], &num_elem_this_blk,
                      &num_nodes_per_elem, &num_attr);
    //std::cout << "elem_blk_id : " << elem_blk_id
    //          << "\nelem_type : " << std::string(elem_type)
    //          << "\nnum_elem_this_blk : " << num_elem_this_blk
    //          << "\nnum_nodes_per_elem : " << num_nodes_per_elem
    //          << "\nnum_attr : " << num_attr
    //          << std::endl;
    _exErr = ex_put_elem_conn(fid, _elmBlock[ieb].id,
                              &(_elmBlock[ieb].conn[0]));
    wrnErrMsg(_exErr, "Problem during writing elementBlock connectivities.");
    // write element attribute
    /*std::vector<double> elmAttrib;
    elmAttrib.resize(_elmBlock[ieb].nElm, 1);
    _exErr = ex_put_elem_attr(fid, 1, &elmAttrib[0]);
    wrnErrMsg(_exErr, "Problem during writing elementBlock attributes.");*/
  }

  // writing names
  std::vector<char *> names;
  for (auto it = ndeSetNames.begin(); it != ndeSetNames.end(); it++)
    names.push_back((char *) (*it).c_str());
  std::cout << names[0] << std::endl;
  _exErr = ex_put_names(fid, EX_NODE_SET, &names[0]);
  wrnErrMsg(_exErr, "Problem writing nodeSet names.");
  names.clear();
  for (auto it = elmBlockNames.begin(); it != elmBlockNames.end(); it++)
    names.push_back((char *) (*it).c_str());
  _exErr = ex_put_names(fid, EX_ELEM_BLOCK, &names[0]);
  wrnErrMsg(_exErr, "Problem writing elementBlock names.");

  // closing the file
  _exErr = ex_close(fid);
  wrnErrMsg(_exErr, "Problem closing the exodusII file.");
}

void exoMesh::exoPopulate(bool updElmLst)
{
  // for node mergin operation we directly update
  // internal database for node coordinates and element
  // connectivities in blocks, so there is no need to
  // rebuild them from nodeSets and elementBlocks
  if (_avoidPopulation) {
    std::cout << "Not populating " << numNdes << "\n";
    return;
  }

  numNdes = 0;
  xCrds.clear();
  yCrds.clear();
  zCrds.clear();
  numElms = 0;
  ndeSetNames.clear();
  elmBlockNames.clear();
  sideSetNames.clear();
  if (updElmLst) {
    for (auto it = _elmBlock.begin(); it != _elmBlock.end(); it++)
      (it->elmIds).clear();
    glbConn.clear();
    for (auto ib = _elmBlock.begin(); ib != _elmBlock.end(); ib++)
      ib->elmIds.clear();
  }

  // gathering node coordinates and updating node sets
  for (auto it1 = _ndeSet.begin(); it1 != _ndeSet.end(); it1++) {
    // sanity check
    //if (it1->usrNdeIds == true)
    //    wrnErrMsg(-666, "User setting for node ids is not supported yet.");
    ndeSetNames.push_back(it1->name);
    // carry out this only if node coordinates are updated
    // because of new node sets
    if (!(it1->usrNdeIds)) {
      for (auto it2 = (it1->crds).begin(); it2 != (it1->crds).end(); it2++) {
        numNdes++;
        (it1->ndeIds).push_back(numNdes);
        xCrds.push_back((*it2)[0]);
        yCrds.push_back((*it2)[1]);
        zCrds.push_back((*it2)[2]);
      }
    }
  }

  // removing empty node sets
  std::vector<ndeSetType> nnss;
  for (auto it1 = _ndeSet.begin(); it1 != _ndeSet.end(); it1++)
    if (it1->nNde > 0)
      nnss.push_back(*it1);
  _ndeSet = nnss;

  // removing empty element blocks
  std::vector<elmBlockType> nebs;
  for (auto it1 = _elmBlock.begin(); it1 != _elmBlock.end(); it1++)
    if (it1->nElm > 0)
      nebs.push_back(*it1);
  _elmBlock = nebs;

  // upadting element sets and block ids
  int blkId = 0;
  for (auto it1 = _elmBlock.begin(); it1 != _elmBlock.end(); it1++) {
    it1->id = ++blkId;
    elmBlockNames.push_back(it1->name);
    // offseting element connectivities
    if ((it1->ndeIdOffset) != 0) {
      for (auto it2 = (it1->conn).begin(); it2 != (it1->conn).end(); it2++)
        (*it2) += it1->ndeIdOffset;
      it1->ndeIdOffset = 0;
    }
    // updating global connectivity
    if (updElmLst) {
      (*it1).elmIds.clear();
      int nn = it1->ndePerElm;
      for (int elmIdx = 0; elmIdx < (*it1).nElm; elmIdx++) {
        std::vector<int> elmConn;
        for (int idx = elmIdx * nn; idx < (elmIdx + 1) * nn; idx++)
          elmConn.push_back(it1->conn[idx]);
        (*it1).elmIds.push_back(numElms + elmIdx);
        glbConn.push_back(elmConn);
      }
    }
    numElms += it1->nElm;
  }

  if (_isVerbose)
    report();
  _isPopulated = true;
}

void exoMesh::report() const
{
  std::cout << " ----- Exodus II Database Report ----- \n";
  std::cout << "Database : " << _ifname << std::endl;
  std::cout << "Nodes: " << numNdes << std::endl;
  std::cout << "Elements: " << numElms << std::endl;
  std::cout << "Node sets: " << getNumberOfNodeSet() << std::endl;
  std::cout << "Element blocks: " << getNumberOfElementBlock() << std::endl;
  std::cout << "Side sets: " << getNumberOfSideSets() << std::endl;
  std::cout << "  Nde      nNde               Name\n";
  std::cout << "  ---      ----              ------\n";
  for (auto ib = _ndeSet.begin(); ib != _ndeSet.end(); ib++)
    std::cout << std::setw(5) << ib->id
              << std::setw(10) << (*ib).nNde
              << std::setw(20) << (*ib).name
              << std::endl;
  std::cout << "  Blk      nElm     eType               Name\n";
  std::cout << "  ---      ----     -----              ------\n";
  for (auto ib = _elmBlock.begin(); ib != _elmBlock.end(); ib++)
    std::cout << std::setw(5) << ib->id
              << std::setw(10) << (*ib).nElm
              << std::setw(10) << EXOMesh::elmTypeStr(ib->eTpe)
              << std::setw(20) << (*ib).name
              << std::endl;
}

void exoMesh::removeByElmIdLst(int blkIdx, const std::vector<int> &idLst)
{
  // range check
  if (blkIdx > getNumberOfElementBlock() || blkIdx < 0)
    wrnErrMsg(-1, "Block Id is out of the range.");

  // the idLst should not eliminate the old block
  // at least one element should remain
  //if (idLst.size() >= _elmBlock[blkIdx].elmIds.size())
  //{
  //    std::cout << idLst.size() << " " << _elmBlock[blkIdx].elmIds.size() << std::endl;
  //    wrnErrMsg(-1, "Can not remove the entire block!.");
  //}

  // now updating the old block
  _isPopulated = false;

  elmBlockType oeb = _elmBlock[blkIdx];
  elmBlockType neb;
  neb.id = oeb.id;
  neb.name = oeb.name;
  neb.ndePerElm = oeb.ndePerElm;
  neb.eTpe = oeb.eTpe;
  neb.ndeIdOffset = 0;
  neb.nElm = 0;

  // elmIds and connectivities
  int oei = -1;
  for (auto it1 = oeb.elmIds.begin(); it1 != oeb.elmIds.end(); it1++) {
    oei++;
    bool stays = true;
    for (auto it2 = idLst.begin(); it2 != idLst.end(); it2++)
      if ((*it2) == (*it1)) {
        stays = false;
        break;
      }

    if (!stays)
      continue;

    neb.nElm++;
    neb.elmIds.push_back(*it1);
    for (int ni = oei * oeb.ndePerElm; ni < (oei + 1) * oeb.ndePerElm; ni++)
      neb.conn.push_back(oeb.conn[ni]);

  }
  std::cout << "Removed " << oeb.elmIds.size() - neb.elmIds.size()
            << " elements from original"
            << std::endl;
  _elmBlock[blkIdx] = neb;
}


void exoMesh::addElmBlkByElmIdLst(const std::string &name,
                                  std::vector<int> &lst)
{
  // preparing database
  if (!_isPopulated);
  exoPopulate(false);

  _isPopulated = false;
  // create element block
  // Assumption here: All elements in the provided list are in the same
  // element block
  elmBlockType neb;
  neb.id = _elmBlock.size() + 1;
  neb.name = name;
  neb.ndeIdOffset = 0;
  //int blkIdx = findElmBlkIdx(lst[0]); // faste but less accurate
  int blkIdx = findElmLstBlkIdx(lst); // slower, more accurate
  if (blkIdx == -1)
    wrnErrMsg(-1, "Elements ids are not registered.");
  //std::cout << "Block Indx = " << blkIdx << std::endl;
  // now adjust the list and remove elements that are
  // not owned by this block
  std::vector<int> owndLst;
  bool allIn;
  owndLst = lstElmInBlck(blkIdx, lst, allIn);
  if (!allIn) {
    std::cout << owndLst.size() << " Elements are in the block.\n";
    std::cout << "Some of the elements in the list are outside the block.\n";
    lst = owndLst;
  }
  neb.elmIds = lst;
  neb.nElm = lst.size();
  neb.eTpe = getBlockElmType(blkIdx);
  neb.ndePerElm = _elmBlock[blkIdx].ndePerElm;
  for (auto eid = lst.begin(); eid != lst.end(); eid++)
    for (int idx = 0; idx < neb.ndePerElm; idx++)
      neb.conn.push_back(glbConn[*eid][idx]);

  // updating the old element block and adding a new one
  removeByElmIdLst(blkIdx, lst);
  addElmBlk(neb);
  exoPopulate(false);
}

void exoMesh::addNdeSetByNdeIdLst(const std::string &name,
                                  const std::vector<int> &idLst)
{
  // preparing database
  if (!_isPopulated);
  exoPopulate(false);

  _isPopulated = false;
  // create and add new node set
  ndeSetType nns;
  nns.id = _ndeSet.size() + 1; // 1-indexed
  nns.nNde = idLst.size();
  nns.name = name;
  nns.usrNdeIds = true;
  nns.ndeIds = idLst;
  addNdeSet(nns);

  // update exodus information
  exoPopulate(false);
}


void exoMesh::snapNdeCrdsZero(double tol)
{
  int nPnt = 0;
  // loop through node sets
  for (auto it1 = _ndeSet.begin(); it1 != _ndeSet.end(); it1++) {
    bool chk;
    for (auto it2 = (it1->crds).begin(); it2 != (it1->crds).end(); it2++) {
      chk = false;
      for (int idx = 0; idx < 3; idx++)
        if (fabs((*it2)[idx]) <= tol) {
          (*it2)[idx] = 0.0;
          chk = true;
        }
      if (chk) nPnt++;
    }
  }
  std::cout << "Number of points snapped : " << nPnt << std::endl;
}

void exoMesh::removeElmBlkByName(std::string blkName)
{
  if (blkName == "")
    return;
  std::vector<elmBlockType> neb;
  for (auto it1 = _elmBlock.begin(); it1 != _elmBlock.end(); it1++) {
    if ((*it1).name.compare(blkName))
      neb.push_back(*it1);
  }
  _elmBlock = neb;
  exoPopulate(true);
}


int exoMesh::findElmBlkIdx(int elmId) const
{
  int blkId = -1;
  int blkIdx = 0;
  for (auto it1 = _elmBlock.begin(); it1 != _elmBlock.end(); it1++) {
    for (auto it2 = (*it1).elmIds.begin(); it2 != (*it1).elmIds.end(); it2++) {
      if ((*it2) == elmId) {
        blkId = it1->id;
        break;
      }
    }
    if (blkId != -1)
      break;
    blkIdx++;
  }
  if (blkId != -1)
    return (blkIdx);
  else
    return (-1);
}

int exoMesh::findElmLstBlkIdx(const std::vector<int> &elmIds) const
{
  int blkIdx = -1;
  int intfCnt = 0;

  if (elmIds.size() == 0)
    return (blkIdx);

  for (int ib = 0; ib < _elmBlock.size(); ib++) {
    std::vector<int> intfLst = {};
    bool allIn = false;
    intfLst = lstElmInBlck(ib, elmIds, allIn);
    if (allIn)
      return (ib);
    else if (intfLst.size() > intfCnt) {
      intfCnt = intfLst.size();
      blkIdx = ib;
    }
  }
  return (blkIdx);
}


std::vector<int> exoMesh::lstElmInBlck(int blkId,
                                       const std::vector<int> &elmIds,
                                       bool allIn) const
{
  std::vector<int> out;
  allIn = true;
  // compare elmIds with list of the elements within block
  std::set<int> current;
  for (auto it = _elmBlock[blkId].elmIds.begin();
       it != _elmBlock[blkId].elmIds.end(); it++)
    current.insert(*it);
  std::pair<std::set<int>::iterator, bool> ret;
  for (auto it = elmIds.begin(); it != elmIds.end(); it++) {
    ret = current.insert(*it);
    if (!ret.second)
      out.push_back(*it);
    else
      allIn = false;
  }
  return (out);
}


void exoMesh::reset()
{
  _isPopulated = false;
  _ndeSet.clear();
  _elmBlock.clear();
  _sideSet.clear();
  fid = 0;
  numNdes = 0;
  numElms = 0;
  xCrds.clear();
  yCrds.clear();
  zCrds.clear();
  glbConn.clear();
  ndeSetNames.clear();
  elmBlockNames.clear();
  sideSetNames.clear();
}


void exoMesh::read(const std::string &ifname)
{
  if (!ifname.empty())
    _ifname = ifname;

  // before reading all internal data base will be reset
  reset();

  int CPU_word_size, IO_word_size;
  float version;
  CPU_word_size = sizeof(float);
  IO_word_size = 0;
  _exErr = 0;

  /* open EXODUS II files */
  fid = ex_open(_ifname.c_str(), EX_READ, &CPU_word_size, &IO_word_size,
                &version);
  wrnErrMsg((fid > 0 ? 0 : -1), "Problem opening file " + _ifname + "\n");

  int numElmBlk;
  int numNdeSet;
  int numSideSet;

  // parameter inquiry from Exodus file
  int num_props;
  int idum;
  float fdum;
  char cdum;

  _exErr = ex_inquire(fid, EX_INQ_API_VERS, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II API version is " << fdum << std::endl;
  _api_v = fdum;

  _exErr = ex_inquire(fid, EX_INQ_DB_VERS, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II Database version is " << fdum << std::endl;
  _dbs_v = fdum;

  _exErr = ex_inquire(fid, EX_INQ_DIM, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of coordinate dimensions is " << num_props << std::endl;
  if (num_props != 3)
    wrnErrMsg(-1, "Only 3D mesh data is supported!\n");

  _exErr = ex_inquire(fid, EX_INQ_NODES, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of points " << num_props << std::endl;
  numNdes = num_props;

  _exErr = ex_inquire(fid, EX_INQ_ELEM, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of elements " << num_props << std::endl;
  numElms = num_props;

  _exErr = ex_inquire(fid, EX_INQ_ELEM_BLK, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numElmBlk = num_props;
  std::cout << "Number of element blocks " << numElmBlk << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_NODE_SETS, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numNdeSet = num_props;
  std::cout << "Number of node sets " << numNdeSet << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_SIDE_SETS, &num_props, &fdum, &cdum);
  wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numSideSet = num_props;
  std::cout << "Number of side sets " << numSideSet << std::endl;

  // nodal coordinates
  std::vector<float> x, y, z;
  xCrds.resize(numNdes, 0);
  yCrds.resize(numNdes, 0);
  zCrds.resize(numNdes, 0);
  _exErr = ex_get_coord(fid, &xCrds[0], &yCrds[0], &zCrds[0]);
  wrnErrMsg(_exErr, "Problem reading nodal coordinates.\n");

  // ids
  std::vector<int> elem_blk_ids;
  elem_blk_ids.resize(numElmBlk, 0);
  _exErr = ex_get_elem_blk_ids(fid, &elem_blk_ids[0]);
  wrnErrMsg(_exErr, "Problem reading element block ids.\n");
  std::vector<int> nde_set_ids;
  nde_set_ids.resize(numNdeSet, 0);
  _exErr = ex_get_node_set_ids(fid, &nde_set_ids[0]);
  wrnErrMsg(_exErr, "Problem reading node set ids.\n");

  // node set and element block names
  std::string blk_name;
  for (int i = 0; i < numElmBlk; i++) {
    int blk_id = elem_blk_ids[i];
    blk_name.resize(MAX_STR_LENGTH);
    _exErr = ex_get_name(fid, EX_ELEM_BLOCK, blk_id, &blk_name[0]);
    wrnErrMsg(_exErr, "Problem reading nodeSet names.");
    std::string::iterator end_pos = std::remove(blk_name.begin(),
                                                blk_name.end(), '\0');
    blk_name.erase(end_pos, blk_name.end());
    //std::cout << blk_id << " " << blk_name << std::endl;
    elmBlockNames.push_back(blk_name);
  }
  blk_name.clear();
  for (int i = 0; i < numNdeSet; i++) {
    blk_name.resize(MAX_STR_LENGTH);
    int blk_id = nde_set_ids[i];
    _exErr = ex_get_name(fid, EX_NODE_SET, blk_id, &blk_name[0]);
    wrnErrMsg(_exErr, "Problem reading nodeSet names.");
    std::string::iterator end_pos = std::remove(blk_name.begin(),
                                                blk_name.end(), '\0');
    blk_name.erase(end_pos, blk_name.end());
    //std::cout << blk_id << " " << blk_name << std::endl;
    ndeSetNames.push_back(blk_name);
  }

  // node sets
  for (int iNS = 0; iNS < numNdeSet; iNS++) {
    ndeSetType ns;
    ns.id = nde_set_ids[iNS];
    ns.name = ndeSetNames[iNS];
    ns.usrNdeIds = true;
    _exErr = ex_get_node_set_param(fid, ns.id, &(ns.nNde), &idum);
    wrnErrMsg(_exErr, "Problem reading node set parameters.\n");
    ns.ndeIds.resize(ns.nNde, 0);
    _exErr = ex_get_node_set(fid, ns.id, &(ns.ndeIds[0]));
    wrnErrMsg(_exErr, "Problem reading node set ids.\n");
    _ndeSet.push_back(ns);
  }

  // element blocks
  std::string elem_type;
  for (int iEB = 0; iEB < numElmBlk; iEB++) {
    elem_type.resize(MAX_STR_LENGTH);
    elmBlockType eb;
    eb.name = elmBlockNames[iEB];
    eb.id = elem_blk_ids[iEB];
    int num_el_in_blk, num_nod_per_el, num_attr;
    _exErr = ex_get_elem_block(fid, eb.id, &elem_type[0], &num_el_in_blk,
                               &num_nod_per_el, &num_attr);
    wrnErrMsg(_exErr, "Problem reading element block parameters.\n");
    std::string::iterator end_pos = std::remove(elem_type.begin(),
                                                elem_type.end(), '\0');
    elem_type.erase(end_pos, elem_type.end());
    eb.nElm = num_el_in_blk;
    eb.ndePerElm = num_nod_per_el;
    eb.eTpe = elmTypeNum(elem_type);
    // read element connectivity
    eb.conn.resize(num_el_in_blk * num_nod_per_el, 0);
    _exErr = ex_get_elem_conn(fid, eb.id, &(eb.conn[0]));
    wrnErrMsg(_exErr, "Problem reading element block connectivities.\n");
    _elmBlock.push_back(eb);
  }
  report();

  // since all data structures are manually populated from the file
  _isPopulated = true;
}

void exoMesh::mergeNodes(double tol)
{
  // build all point kdtree
  int nDim = 3;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(numNdes, nDim);
  for (int iPnt = 0; iPnt < numNdes; iPnt++) {
    pntCrd[iPnt][0] = xCrds[iPnt];
    pntCrd[iPnt][1] = yCrds[iPnt];
    pntCrd[iPnt][2] = zCrds[iPnt];
  }
  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, numNdes, nDim);

  // finding duplicate node ids
  // exomesh is one-indexed
  std::map<int, int> dupNdeMap;
  int nNib = 2;
  ANNpoint qryPnt;
  ANNidxArray nnIdx = new ANNidx[nNib];
  ANNdistArray dists = new ANNdist[nNib];
  qryPnt = annAllocPt(nDim);
  for (int iNde = 0; iNde < numNdes; iNde++) {
    qryPnt[0] = xCrds[iNde];
    qryPnt[1] = yCrds[iNde];
    qryPnt[2] = zCrds[iNde];
    kdTree->annkSearch(qryPnt, nNib, nnIdx, dists);
    if (dists[1] <= tol)
      dupNdeMap[nnIdx[0] + 1] = nnIdx[1] + 1;
  }
  std::cout << "Found " << dupNdeMap.size() << " duplicate nodes.\n";

  // removing duplicated nodal coordinates and
  // creating map from old ids to new ones
  std::vector<double> xn, yn, zn;
  std::map<int, int> old2NewNde;
  int newNdeIdx = 0;
  for (int iNde = 0; iNde < numNdes; iNde++) {
    auto ite = dupNdeMap.find(iNde + 1);
    if (ite == dupNdeMap.end()) {
      xn.push_back(xCrds[iNde]);
      yn.push_back(yCrds[iNde]);
      zn.push_back(zCrds[iNde]);
    }
  }
  // mapping
  delete kdTree;
  for (int iPnt = 0; iPnt < xn.size(); iPnt++) {
    pntCrd[iPnt][0] = xn[iPnt];
    pntCrd[iPnt][1] = yn[iPnt];
    pntCrd[iPnt][2] = zn[iPnt];
  }
  kdTree = new ANNkd_tree(pntCrd, xn.size(), nDim);
  for (int iNde = 0; iNde < numNdes; iNde++) {
    qryPnt[0] = xCrds[iNde];
    qryPnt[1] = yCrds[iNde];
    qryPnt[2] = zCrds[iNde];
    kdTree->annkSearch(qryPnt, 1, nnIdx, dists);
    if (dists[0] <= tol)
      old2NewNde[iNde + 1] = nnIdx[0] + 1;
    else {
      cerr << "Found a node in database which is not copied properly \n";
      exit(1);
    }

  }

  // sanity check and node coordinate re-assignment
  int max, min;
  max = 0;
  min = INT_MAX;
  for (auto im = old2NewNde.begin(); im != old2NewNde.end(); im++) {
    max = std::max(max, im->second);
    min = std::min(min, im->second);
  }
  std::cout << "Max/Min : " << max << "/" << min << std::endl;

  xCrds.clear();
  yCrds.clear();
  zCrds.clear();
  xCrds.assign(xn.begin(), xn.end());
  yCrds.assign(yn.begin(), yn.end());
  zCrds.assign(zn.begin(), zn.end());
  cout << "Number of nodes changed from " << numNdes
       << " to " << xn.size() << "\n";
  numNdes = xn.size();

  // update nodesets node ids
  // update nodeset node coordinates
  std::set<int> idTally;
  std::set<int> idTally2;
  for (auto ins = _ndeSet.begin(); ins != _ndeSet.end(); ins++) {

    std::set<int> idsToKeep;
    int nsLocalNdeIds = 0;
    for (auto itc = (ins->crds).begin(); itc != (ins->crds).end(); itc++) {
      qryPnt[0] = (*itc)[0];
      qryPnt[1] = (*itc)[1];
      qryPnt[2] = (*itc)[2];
      kdTree->annkSearch(qryPnt, 1, nnIdx, dists);
      auto ret = idTally.insert(nnIdx[0]);
      if (ret.second)
        idsToKeep.insert(nsLocalNdeIds++);
    }
    cout << "Original node set with " << ins->nNde
         << " nodes after cleaning becomes " << idsToKeep.size() << "\n";

    std::vector<std::vector<double> > ncrds;
    std::vector<double> crds;
    crds.resize(3, -1);
    for (auto k = idsToKeep.begin(); k != idsToKeep.end(); k++) {
      crds[0] = ((ins->crds)[*k])[0];
      crds[1] = ((ins->crds)[*k])[1];
      crds[2] = ((ins->crds)[*k])[2];
      ncrds.push_back(crds);
    }
    cout << "Size of node set coords " << ncrds.size() << "\n";

    // updating node indices
    std::set<int> newNSNdeIds;
    for (auto itn = (ins->ndeIds).begin(); itn != (ins->ndeIds).end(); itn++) {
      int nid = *itn;
      int newId = old2NewNde[nid];
      auto ret = idTally2.insert(newId);
      if (ret.second)
        newNSNdeIds.insert(newId);
    }
    (ins->ndeIds).clear();
    (ins->ndeIds).insert((ins->ndeIds).end(),
                         newNSNdeIds.begin(), newNSNdeIds.end());
    cout << "NdeId size after " << (ins->ndeIds).size() << "\n";

    // updating the nodeset data structure
    (ins->nNde) = idsToKeep.size();
    (ins->crds) = ncrds;

  }


  // update element block connectivities
  max = 0;
  min = INT_MAX;
  for (auto ieb = _elmBlock.begin(); ieb != _elmBlock.end(); ieb++) {
    for (auto itn = (ieb->conn).begin(); itn != (ieb->conn).end(); itn++) {
      int nid = *itn;
      *itn = old2NewNde[nid];
      min = std::min(min, *itn);
      max = std::max(max, *itn);
    }
  }
  std::cout << "Min conn = " << min << std::endl;
  std::cout << "Max conn = " << max << std::endl;

  // sidesets does not need to be updated since
  // they contain only element ids

  // TODO: avoid using this flag
  // avoiding database rebuilding since it is already upto
  // date
  _avoidPopulation = true;

  // deleting KDTree
  if (kdTree)
    delete kdTree;
  delete nnIdx;
  delete dists;
}

#endif
