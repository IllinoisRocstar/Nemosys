#include "meshBase.H"
#include "exoMesh.H"
#include <fstream>
#include <algorithm>
#include "exodusII.h"

using namespace EXOMesh;

VTKCellType EXOMesh::e2vEMap ( elementType et)
{
  std::map<elementType,VTKCellType> eMap = 
  {
    {elementType::QUAD, VTK_QUAD},
    {elementType::TRIANGLE, VTK_TRIANGLE},
    {elementType::TETRA, VTK_TETRA},
    {elementType::WEDGE, VTK_WEDGE},
    {elementType::HEX, VTK_HEXAHEDRON},
    {elementType::OTHER, VTK_EMPTY_CELL}
  };
  return eMap[et];
}

elementType EXOMesh::v2eEMap (VTKCellType vt)
{
  std::map<VTKCellType,elementType> eMap = 
  {
    {VTK_QUAD, elementType::QUAD},
    {VTK_QUADRATIC_QUAD, elementType::QUAD},
    {VTK_TRIANGLE, elementType::TRIANGLE},
    {VTK_QUADRATIC_TRIANGLE, elementType::TRIANGLE},
    {VTK_HEXAHEDRON, elementType::HEX},
    {VTK_QUADRATIC_HEXAHEDRON, elementType::HEX},
    {VTK_TETRA, elementType::TETRA},
    {VTK_QUADRATIC_TETRA, elementType::TETRA},
    {VTK_WEDGE, elementType::WEDGE},
    {VTK_EMPTY_CELL, elementType::OTHER}
  };
  return eMap[vt];
}

surfaceBCTag EXOMesh::bcTagNum(std::string& tag)
{
  std::transform(tag.begin(), tag.end(), tag.begin(), ::tolower);
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

elementType EXOMesh::elmTypeNum (std::string tag) 
{
  std::transform(tag.begin(), tag.end(), tag.begin(), ::tolower);
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

std::string EXOMesh::elmTypeStr (elementType tag)
{
  if (tag == elementType::QUAD) return "QUAD";
  if (tag == elementType::TRIANGLE) return "TRIANGLE";
  if (tag == elementType::HEX) return "HEX";
  if (tag == elementType::TETRA) return "TETRA";
  if (tag == elementType::WEDGE) return "WEDGE";
  return "OTHER";
}

int EXOMesh::elmNumNde (elementType tag, int order)
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

int EXOMesh::elmNumSrf (elementType tag)
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
    _isOpen(true)
{}

exoMesh::~exoMesh() 
{
    if (_isOpen)
        ex_close(fid);
}
    
    
void EXOMesh::wrnErrMsg(int errCode, std::string msg)
{
    if (errCode<0) {
        std::cerr<< "Error: " << msg 
                 << " with error code "
                 << errCode << std::endl;
        throw;
    } 
    else if (errCode>0)
    {
        std::cerr<< "Warning: code number " << errCode << std::endl;
    }
}

void exoMesh::write()
{
    // preparing database
    // regadless we update it
    exoPopulate(true);

    // writing to file
    int comp_ws = 8;
    int io_ws = 8;
    fid = ex_create(_ifname.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

    // initializing exodus database
    _exErr = ex_put_init(fid,"Nemosys ExodusII database", 3, numNdes, numElms, 
            _elmBlock.size(), _ndeSet.size(), _sideSet.size() );
    wrnErrMsg(_exErr, "Problem during initialization of exodusII database");

    // writing node coordinates
    _exErr = ex_put_coord(fid, &xCrds[0], &yCrds[0], &zCrds[0]); 
    wrnErrMsg(_exErr, "Problem during writing node coordinates.");
    
    // writing node sets
    for (int ins=0; ins<_ndeSet.size(); ins++)
    {
        _exErr = ex_put_node_set_param(fid, _ndeSet[ins].id, _ndeSet[ins].nNde, 0);
        wrnErrMsg(_exErr, "Problem during writing nodeSet parameteres.");
        _exErr = ex_put_node_set(fid, _ndeSet[ins].id, &(_ndeSet[ins].ndeIds[0]) );
        wrnErrMsg(_exErr, "Problem during writing nodeSet node Ids.");
    }

    // writing element blocks
    for (int ieb=0; ieb<_elmBlock.size(); ieb++)
    {
        std::string eTpe = elmTypeStr(_elmBlock[ieb].eTpe);
        _exErr = ex_put_elem_block(fid, _elmBlock[ieb].id, eTpe.c_str(), 
                _elmBlock[ieb].nElm, _elmBlock[ieb].ndePerElm, 1);
        wrnErrMsg(_exErr, "Problem during writing elementBlock parameteres.");
        int elem_blk_id, num_elem_this_blk, num_nodes_per_elem, num_attr;
        char elem_type[MAX_STR_LENGTH+1];
        elem_blk_id = ieb+1;
        ex_get_elem_block (fid, elem_blk_id, elem_type, &num_elem_this_blk,
            &num_nodes_per_elem, &num_attr);
        //std::cout << "elem_blk_id : " << elem_blk_id 
        //          << "\nelem_type : " << std::string(elem_type)
        //          << "\nnum_elem_this_blk : " << num_elem_this_blk
        //          << "\nnum_nodes_per_elem : " << num_nodes_per_elem
        //          << "\nnum_attr : " << num_attr
        //          << std::endl;
        //end testing
        _exErr = ex_put_elem_conn(fid, _elmBlock[ieb].id, &(_elmBlock[ieb].conn[0]));
        wrnErrMsg(_exErr, "Problem during writing elementBlock connectivities.");
        // write element attribute
        std::vector<double> elmAttrib;
        elmAttrib.resize(_elmBlock[ieb].nElm, 1);
        _exErr = ex_put_elem_attr(fid, 1, &elmAttrib[0]);
        wrnErrMsg(_exErr, "Problem during writing elementBlock attributes.");
    }

    // writing names
    _exErr = ex_put_names(fid, EX_NODE_SET, &ndeSetNames[0]);
    wrnErrMsg(_exErr, "Problem writing nodeSet names.");
    _exErr = ex_put_names(fid, EX_ELEM_BLOCK, &elmBlockNames[0]);
    wrnErrMsg(_exErr, "Problem writing elementBlock names.");

    // closing the file
    _exErr = ex_close(fid);
    wrnErrMsg(_exErr, "Problem closing the exodusII file.");
}

void exoMesh::exoPopulate(bool updElmLst) 
{
    numNdes = 0;
    numElms = 0;
    xCrds.clear();
    yCrds.clear();
    zCrds.clear();
    ndeSetNames.clear();
    elmBlockNames.clear();
    sideSetNames.clear();
    if (updElmLst)
    {
        glbConn.clear();
        for (auto ib=_elmBlock.begin(); ib!=_elmBlock.end(); ib++)
            ib->elmIds.clear();
    }

    // gathering node coordinates and updating node sets
    for (auto it1=_ndeSet.begin(); it1!=_ndeSet.end(); it1++)
    {
        // sanity check
        if (it1->usrNdeIds == true)
            wrnErrMsg(-666, "User setting for node ids is not supported yet.");
        ndeSetNames.push_back(const_cast<char*>((it1->name).c_str()));
        for (auto it2=(it1->crds).begin(); it2!=(it1->crds).end(); it2++)
        {
            numNdes++;
            (it1->ndeIds).push_back(numNdes);
            xCrds.push_back((*it2)[0]);
            yCrds.push_back((*it2)[1]);
            zCrds.push_back((*it2)[2]);
        }
    }

    // removing empty element blocks
    std::vector<elmBlockType> nebs;
    for (auto it1=_elmBlock.begin(); it1!=_elmBlock.end(); it1++)
        if (it1->nElm > 0)
            nebs.push_back(*it1);
    _elmBlock = nebs;

    // upadting element sets and block ids
    int blkId = 0;
    for (auto it1=_elmBlock.begin(); it1!=_elmBlock.end(); it1++)
    {
        it1->id = ++blkId;
        elmBlockNames.push_back(const_cast<char*>((it1->name).c_str()));        
        // offseting element connectivities
        if ( (it1->ndeIdOffset) != 0 )
        {
            std::cout << "Non-zero node offset\n";
            for (auto it2=(it1->conn).begin(); it2!=(it1->conn).end(); it2++)
                (*it2) += it1->ndeIdOffset;
            it1->ndeIdOffset = 0;
        }
        // updating global connectivity
        if (updElmLst)
        {
            int nn = it1->ndePerElm;
            for (int elmIdx=0; elmIdx<(*it1).nElm; elmIdx++)
            {
                std::vector<int> elmConn;
                for (int idx=elmIdx*nn; idx<(elmIdx+1)*nn; idx++)
                    elmConn.push_back(it1->conn[idx]);       
                (*it1).elmIds.push_back(numElms+elmIdx);
                glbConn.push_back(elmConn);
            }
        }
        numElms += it1->nElm;
    }

    report();
    _isPopulated = true;
}

void exoMesh::report() const
{
    std::cout << " ----- Exodus II Database Statistics ----- \n";
    std::cout << "#Nodes = " << numNdes << std::endl;
    std::cout << "#Elements = " << numElms << std::endl;
    std::cout << "#Node sets = " << getNumberOfNodeSet() << std::endl;
    std::cout << "#Element blocks = " << getNumberOfElementBlock() << std::endl;
    std::cout << "#Side sets = " << getNumberOfSideSets() << std::endl;
    std::cout << "  Blk      nElm     eType\n";
    std::cout << "  ---      ----     -----\n";
    for (auto ib=_elmBlock.begin(); ib!=_elmBlock.end(); ib++)
        std::cout << std::setw(5) << ib->id
                  << std::setw(10) << (*ib).elmIds.size()
                  << std::setw(10) << EXOMesh::elmTypeStr(ib->eTpe)
                  << std::endl;
}

void exoMesh::removeByElmIdLst(int blkIdx, std::vector<int>& idLst)
{
    // range check
    if (blkIdx>getNumberOfElementBlock() || blkIdx<0) 
        wrnErrMsg(-1, "Block Id is out of the range.");

    // the idLst should not eliminate the old block
    // at least one element should remain
    //if (idLst.size() >= _elmBlock[blkIdx].elmIds.size())
    //{
    //    std::cout << idLst.size() << " " << _elmBlock[blkIdx].elmIds.size() << std::endl; 
    //    wrnErrMsg(-1, "Can not remove the entire block!.");
    //}

    // now updating the old block
    _isPopulated=false;

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
    for (auto it1=oeb.elmIds.begin(); it1!=oeb.elmIds.end(); it1++)
    {
        oei++;
        bool stays = true;
        for (auto it2=idLst.begin(); it2!=idLst.end(); it2++)
            if ( (*it2) == (*it1) )
            {
                stays = false;
                break;
            }

        if (!stays)
            continue;

        neb.nElm++;
        neb.elmIds.push_back(*it1);
        for (int ni=oei*oeb.ndePerElm; ni<(oei+1)*oeb.ndePerElm; ni++)
            neb.conn.push_back(oeb.conn[ni]);

    }
    _elmBlock[blkIdx] = neb;
}


void exoMesh::addElmBlkByElmIdLst(std::string name, std::vector<int>& lst)
{ 
    // preparing database
    if (!_isPopulated);
        exoPopulate(false);

    _isPopulated=false;

    // create element block
    // Assumption here: All elements in the provided list are in the same 
    // element block
    elmBlockType neb;
    neb.id = _elmBlock.size()+1; 
    neb.nElm = lst.size();
    neb.name = name;
    neb.ndeIdOffset = 0;
    neb.elmIds = lst;
    int blkIdx = findElmBlkIdx(lst[0]);
    if (blkIdx == -1) 
        wrnErrMsg(-1, "Elements ids are not registered.");
    //std::cout << "Block Indx = " << blkIdx << std::endl;
    neb.eTpe = getBlockElmType(blkIdx); 
    neb.ndePerElm = _elmBlock[blkIdx].ndePerElm;
    for (auto eid=lst.begin(); eid!=lst.end(); eid++)
        for (int idx=0; idx<neb.ndePerElm; idx++)
            neb.conn.push_back( glbConn[*eid][idx] );
    
    // updating the old element block and adding a new one
    removeByElmIdLst(blkIdx, lst);
    addElmBlk(neb);
    exoPopulate(false);
}

int exoMesh::findElmBlkIdx(int elmId) const
{
    int blkId = -1;
    int blkIdx = 0;
    for (auto it1=_elmBlock.begin(); it1!=_elmBlock.end(); it1++)
    {
        for (auto it2=(*it1).elmIds.begin(); it2!=(*it1).elmIds.end(); it2++)
        {
            if ( (*it2)==elmId )
            {
                blkId = it1->id;
                break;
            }
        }
        if ( blkId != -1)
            break;
        blkIdx++;
    }
    if (blkId != -1)
        return(blkIdx);
    else
        return(-1);
}
