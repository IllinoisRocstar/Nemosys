#include "meshBase.H"
#include "pntMesh.H"
#include <fstream>
#include <algorithm>


using namespace PNTMesh;

VTKCellType PNTMesh::p2vEMap ( elementType et)
{
  std::map<elementType,VTKCellType> eMap = 
  {
    {elementType::BAR, VTK_LINE},
    {elementType::QUADRILATERAL, VTK_QUAD},
    {elementType::TRIANGLE, VTK_TRIANGLE},
    {elementType::HEXAGON, VTK_HEXAHEDRON},
    {elementType::SPHERICAL, VTK_EMPTY_CELL},
    {elementType::CYLINDRICAL, VTK_EMPTY_CELL},
    {elementType::BRICK, VTK_HEXAHEDRON},
    {elementType::LAGRANGE_BRICK, VTK_HEXAHEDRON},
    {elementType::TETRAHEDRON, VTK_TETRA},
    {elementType::HEXPRISM, VTK_HEXAGONAL_PRISM},
    {elementType::PRISMATIC, VTK_WEDGE},
    {elementType::OTHER, VTK_EMPTY_CELL}
  };
  return eMap[et];
}

elementType PNTMesh::v2pEMap (VTKCellType vt)
{
  std::map<VTKCellType,elementType> eMap = 
  {
    {VTK_LINE, elementType::BAR},
    {VTK_QUAD, elementType::QUADRILATERAL},
    {VTK_QUADRATIC_QUAD, elementType::QUADRILATERAL},
    {VTK_TRIANGLE, elementType::TRIANGLE},
    {VTK_QUADRATIC_TRIANGLE, elementType::TRIANGLE},
    {VTK_HEXAHEDRON, elementType::HEXAGON},
    {VTK_QUADRATIC_HEXAHEDRON, elementType::HEXAGON},
    {VTK_TETRA, elementType::TETRAHEDRON},
    {VTK_QUADRATIC_TETRA, elementType::TETRAHEDRON},
    {VTK_HEXAGONAL_PRISM, elementType::HEXPRISM},
    {VTK_WEDGE, elementType::PRISMATIC},
    {VTK_EMPTY_CELL, elementType::OTHER}
  };
  return eMap[vt];
}

surfaceBCTag PNTMesh::bcTagNum(std::string& tag)
{
  std::transform(tag.begin(), tag.end(), tag.begin(), ::tolower);
  if (tag == "reflective") return surfaceBCTag::REFLECTIVE;
  if (tag == "void") return surfaceBCTag::VOID;
  std::cerr << "Unknown surface tag " << tag << std::endl;
  throw;
}

std::string PNTMesh::bcTagStr(int tag) 
{
  if (tag == REFLECTIVE) return "REFLECTIVE";
  if (tag == VOID) return "VOID";
  std::cerr << "Unknown surface tag " << tag << std::endl;
  throw;
}

elementType PNTMesh::elmTypeNum (std::string tag) 
{
  std::transform(tag.begin(), tag.end(), tag.begin(), ::tolower);
  if (tag == "bar") return elementType::BAR;
  if (tag == "quadrilateral") return elementType::QUADRILATERAL;
  if (tag == "triangle") return elementType::TRIANGLE;
  if (tag == "hexagon") return elementType::HEXAGON;
  if (tag == "spherical") return elementType::SPHERICAL;
  if (tag == "cylindrical") return elementType::CYLINDRICAL;
  if (tag == "brick") return elementType::BRICK;
  if (tag == "lagrange_brick") return elementType::LAGRANGE_BRICK;
  if (tag == "tetrahedron") return elementType::TETRAHEDRON;
  if (tag == "hexprism") return elementType::HEXPRISM;
  if (tag == "prismatic") return elementType::PRISMATIC;
  std::cout 
      << "Warning : Element type "
      << tag
      << " may be not supported\n";
  return elementType::OTHER;
}

std::string PNTMesh::elmTypeStr (elementType tag)
{
  if (tag == elementType::BAR) return "BAR";
  if (tag == elementType::QUADRILATERAL) return "QUADRILATERAL";
  if (tag == elementType::TRIANGLE) return "TRIANGLE";
  if (tag == elementType::HEXAGON) return "HEXAGON";
  if (tag == elementType::SPHERICAL) return "SPHERICAL";
  if (tag == elementType::CYLINDRICAL) return "CYLINDRICAL";
  if (tag == elementType::BRICK) return "BRICK";
  if (tag == elementType::LAGRANGE_BRICK) return "LAGRANGE_BRICK";
  if (tag == elementType::TETRAHEDRON) return "TETRAHEDRON";
  if (tag == elementType::HEXPRISM) return "HEXPRISM";
  if (tag == elementType::PRISMATIC) return "PRISMATIC";
  return "OTHER";
}

int PNTMesh::elmNumNde (elementType tag, int order)
{
  if (tag == elementType::BAR) return (2+order-1);
  if (tag == elementType::QUADRILATERAL && order == 1) return 4;
  if (tag == elementType::QUADRILATERAL && order == 2) return 8;
  if (tag == elementType::TRIANGLE && order == 1) return 3;
  if (tag == elementType::TRIANGLE && order == 2) return 6;
  if (tag == elementType::HEXAGON && order == 1) return 8;
  if (tag == elementType::HEXAGON && order == 2) return 27;
  if (tag == elementType::BRICK && order == 1) return 8;
  if (tag == elementType::BRICK && order == 2) return 27;
  if (tag == elementType::LAGRANGE_BRICK && order == 1) return 8;
  if (tag == elementType::LAGRANGE_BRICK && order == 2) return 27;
  if (tag == elementType::TETRAHEDRON && order == 1) return 4;
  if (tag == elementType::TETRAHEDRON && order == 2) return 10;
  if (tag == elementType::HEXPRISM && order == 1) return 16;
  if (tag == elementType::PRISMATIC && order == 1) return 6;
  if (tag == elementType::PRISMATIC && order == 2) return 15;
  std::cerr << "Unknown element/order combination\n";
  throw;
}

int PNTMesh::elmNumSrf (elementType tag)
{
  if (tag == elementType::BAR) return 2;
  if (tag == elementType::QUADRILATERAL) return 4;
  if (tag == elementType::TRIANGLE) return 3;
  if (tag == elementType::HEXAGON) return 8;
  if (tag == elementType::BRICK) return 6;
  if (tag == elementType::LAGRANGE_BRICK) return 6;
  if (tag == elementType::TETRAHEDRON) return 4;
  if (tag == elementType::HEXPRISM) return 10;
  if (tag == elementType::PRISMATIC) return 6;
  std::cerr << "Unknown element type\n";
  throw;
}
//////////////////////////////////
// pntMesh class 
//////////////////////////////////

pntMesh::pntMesh() : 
    ifname(""), isSupported(true)
{}

pntMesh::pntMesh(std::string ifname) :
    ifname(ifname), isSupported(true)
{
    
  // NOTE: If file is very large different means are needed
  // current implementation is focused on the speed of reading
  // information from the file.
  
  std::ifstream fs(ifname);
  if (!fs.good())
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }

  // getting size
  fs.seekg(0,std::ios::end);
  std::streampos length = fs.tellg();
  fs.seekg(0,std::ios::beg);
    
  // read entire file into buffer and then pass to stream
  std::vector<char>  buff(length);
  fs.read(&buff[0],length);
  std::stringstream ss;
  ss.rdbuf()->pubsetbuf(&buff[0],length);

  /////////////////////////////
  // processing CARD 01
  // //////////////////////////
  // header
  int dmi;
  std::string dms;

  ss >> numVertices;
  ss >> numElements; 
  ss >> numDimensions;
  ss >> numBlocks; 
  ss >> numSurfaces; 
  ss >> numSurfInternal; 
  ss >> numSurfBoundary;
  ss >> dmi; // dummy
  ss >> dmi;
  ss >> dmi;

  std::cout << "Reading PNT Mesh..." << std::endl;
  std::cout << "Number of vertices : " << numVertices << std::endl;
  std::cout << "Number of eLements : " << numElements << std::endl;
  std::cout << "Number of dimensions : " << numDimensions << std::endl;
  std::cout << "Number of blocks : " << numBlocks << std::endl;
  std::cout << "Number of surfaces : " << numSurfaces << std::endl;
  std::cout << "Number of internals surface : " << numSurfInternal << std::endl;
  std::cout << "Number of boundary surface  : " << numSurfBoundary << std::endl;

  /////////////////////////////
  // processing CARD 02
  // //////////////////////////
  // nodal coords
  pntCrds.resize(numVertices, std::vector<double>(3,0));
  for (int id=0; id<numDimensions; id++)
  {
    for (int iv=0; iv<numVertices; iv++)
    {
      ss >> pntCrds[iv][id];
    }
  }

  /////////////////////////////
  // processing CARD 03
  // //////////////////////////
  // block data
  elmBlks.resize(numBlocks); 
  for (int iBlk=0; iBlk<numBlocks; iBlk++)
  {
    blockType nb;

    // block header
    ss >> nb.numElementsInBlock;
    ss >> nb.numBoundarySurfacesInBlock;
    ss >> nb.nodesPerElement;
    ss >> dmi;
    ss >> nb.ordIntrp;
    ss >> nb.ordEquat;
    ss >> dms; std::cout<<dms<<std::endl; nb.eTpe = elmTypeNum(dms);
    ss >> nb.regionName;

    std::cout << "================================================" << std::endl;
    std::cout << "Processing block " << nb.regionName << std::endl;
    std::cout << "Num Elements : " << nb.numElementsInBlock << std::endl;
    std::cout << "Num Boundary Surfaces : " << nb.numBoundarySurfacesInBlock << std::endl;
    std::cout << "Num Nodes/Element : " << nb.nodesPerElement << std::endl;
    std::cout << "Element order : " << nb.ordIntrp << std::endl;
    std::cout << "Equation order : " << nb.ordEquat << std::endl;
    std::cout << "Element Type : " << elmTypeStr(nb.eTpe) << std::endl;

    // make sure all blocks are tetrahedral or triangular
    if (nb.eTpe != elementType::TETRAHEDRON && 
        nb.eTpe != elementType::TRIANGLE) 
        isSupported=false;

    // element connectivity array
    std::cout << "Reading block...";
    nb.eConn.resize(nb.numElementsInBlock);
    for (int ibe=0; ibe<nb.numElementsInBlock; ibe++)
    {
      std::vector<idTyp> conn;
      conn.resize(nb.nodesPerElement);
      for (int ine=0; ine<nb.nodesPerElement; ine++)
      {
        ss >> conn[ine];
        --conn[ine]; // pnt mesh is 1-indexed
      }
      nb.eConn[ibe] = conn;
      elmConn.push_back(conn);
      elmTyp.push_back(nb.eTpe);
      elmOrd.push_back(nb.ordIntrp);
    }

    // surface boundary tags
    std::cout << "...";
    nb.srfBCTag.resize(nb.numBoundarySurfacesInBlock);
    for (int ist=0; ist<nb.numBoundarySurfacesInBlock; ist++)
    {
      ss >> dms; 
      nb.srfBCTag[ist] = bcTagNum(dms);
    }
   
    // element surface reference number
    std::cout << "...";
    nb.srfBCEleRef.resize(2*nb.numBoundarySurfacesInBlock);
    for (int ise=0; ise<nb.numBoundarySurfacesInBlock; ise++)
    {
      idTyp er;
      ss >> er;
      nb.srfBCEleRef.push_back(er);
      ss >> er;
      nb.srfBCEleRef.push_back(er);
    }

    // adjacancy header
    std::cout << "...";
    ss >> dmi;
    ss >> nb.numSurfPerEleInBlock;
    int ni = nb.numSurfPerEleInBlock*nb.numElementsInBlock;
    nb.glbSrfId.resize(ni);
    nb.adjBlkId.resize(ni);
    nb.adjElmId.resize(ni);
    nb.adjRefId.resize(ni);

    // global surface ids 
    std::cout << "...";
    for (int i=0; i<ni; i++)
      ss >> nb.glbSrfId[i];
    
    // adjacent block ids
    std::cout << "...";
    for (int i=0; i<ni; i++)
      ss >> nb.adjBlkId[i];
    
    // adjacent element ids
    std::cout << "...";
    for (int i=0; i<ni; i++)
      ss >> nb.adjElmId[i];
    
    // adjacent reference ids
    std::cout << "...\n";
    for (int i=0; i<ni; i++)
      ss >> nb.adjRefId[i];

    // block finishes
    elmBlks[iBlk] = nb; 
  
  }
 
  // closing the file
  fs.close();
}

// constructor from meshBase object
// reads a meshbase object and populates internal database
// assumes all elements of certain type reside in separate
// blocks,
// TODO: Should be exteneded to get block data, BCs, etc.
// and prepare internal datastructure properly based on these
// information
pntMesh::pntMesh(const meshBase* imb, int dim, int nBlk, BlockMap& elmBlkMap) :
    numDimensions(dim), numBlocks(nBlk), numSurfaces(0), numSurfInternal(0), numSurfBoundary(0)
{
  numVertices = imb->getNumberOfPoints();
  numElements = imb->getNumberOfCells();

  // populating point coordinates
  pntCrds.resize(numVertices);
  for (int ipt=0; ipt<numVertices; ipt++)
  {
    std::vector<double> crds;
    crds.resize(3,0.);
    imb->getDataSet()->GetPoint(ipt, &crds[0]);
    pntCrds[ipt] = crds;
  }

  // populating cell connectivity
  elmConn.resize(numElements);
  elmTyp.resize(numElements);
  for (int i = 0; i < numElements; ++i)
  {
    vtkIdList* point_ids = imb->getDataSet()->GetCell(i)->GetPointIds();
    int numComponent = point_ids->GetNumberOfIds();
    std::vector<int> cn;
    cn.resize(numComponent);
    for (int j = 0; j < numComponent; ++j)
      cn[j] = point_ids->GetId(j);
    elmConn[i] = cn;
    
    VTKCellType type_id = static_cast<VTKCellType>(imb->getDataSet()->GetCellType(i));
    elmTyp[i] = v2pEMap(type_id);
  }

  // populating element block ids, and order
  elmBlkId.resize(numElements,-1);
  elmLocalId.resize(numElements,-1);
  elmOrd.resize(numElements, -1);
  for (int iBlk=0; iBlk<numBlocks; iBlk++)
  {
    int lid = 0;
    for (auto it = elmBlkMap[iBlk].elmIds.begin();
              it != elmBlkMap[iBlk].elmIds.end();
              it++)
    {
        elmBlkId[*it] = iBlk;
        elmLocalId[*it] = lid++;
        elmOrd[*it] = elmBlkMap[iBlk].ordIntrp;
    }
  }

  // analyzing element blocks and populating some of the
  // fileds
  elmBlks = elmBlkMap;
  for (int iBlk=0; iBlk<numBlocks; iBlk++)
  {
    elmBlks[iBlk].numElementsInBlock = elmBlks[iBlk].elmIds.size();
    elmBlks[iBlk].nodesPerElement = elmNumNde(elmBlks[iBlk].eTpe, elmBlks[iBlk].ordIntrp);
    elmBlks[iBlk].numSurfPerEleInBlock = elmNumSrf(elmBlks[iBlk].eTpe);
    std::cout << "Block Name = " << elmBlks[iBlk].regionName << std::endl;
    std::cout << "Block size = " << elmBlks[iBlk].numElementsInBlock << std::endl;
  }

  //populate remaining fields
  pntPopulate(imb);

}

int pntMesh::getNumberOfPoints() const
{
  return numVertices;
}

int pntMesh::getNumberOfCells() const
{
  return numElements; 
}

std::vector<double> pntMesh::getPointCrd(int id) const
{
  return pntCrds[id];
}

std::vector<int> pntMesh::getElmConn(int id) const
{
  return elmConn[id];
}

std::vector<int> pntMesh::getElmConn(int id, VTKCellType vct) const
{
    std::vector<int> ci = getElmConn(id);
    std::vector<int> co;
    co = ci;
    if (vct==VTK_QUADRATIC_TRIANGLE)
    {
        co[0] = ci[0];
        co[1] = ci[2];
        co[2] = ci[4];
        co[3] = ci[1];
        co[4] = ci[3];
        co[5] = ci[5];
    }
    else if (vct==VTK_QUADRATIC_HEXAHEDRON)
    {
        co[0] = ci[0];
        co[1] = ci[2];
        co[2] = ci[8];
        co[3] = ci[6];
        co[4] = ci[18];
        co[5] = ci[20];
        co[6] = ci[26];
        co[7] = ci[24];
        co[8] = ci[1];
        co[9] = ci[5];
        co[10] = ci[7];
        co[11] = ci[3];
        co[12] = ci[19];
        co[13] = ci[23];
        co[14] = ci[25];
        co[15] = ci[21];
        co[16] = ci[9];
        co[17] = ci[11];
        co[18] = ci[17];
        co[19] = ci[15];
        co[20] = ci[12];
        co[21] = ci[14];
        co[22] = ci[10];
        co[23] = ci[16];
        co[24] = ci[4];
        co[25] = ci[22];
        co[26] = ci[13];
    }
    return co;
}

std::vector<int> pntMesh::getPntConn(std::vector<int>& ci, elementType et, int eo) const
{
    std::vector<int> co;
    co = ci;
    if (et == TRIANGLE && eo == 2)
    {
        co[0] = ci[0];
        co[2] = ci[1];
        co[4] = ci[2];
        co[1] = ci[3];
        co[3] = ci[4];
        co[5] = ci[5];
    }
    else if ((et == HEXAGON || et == LAGRANGE_BRICK || et == BRICK) 
            && eo == 2)
    {
        co[0] = ci[0];
        co[2] = ci[1];
        co[8] = ci[2];
        co[6] = ci[3];
        co[18] = ci[4];
        co[20] = ci[5];
        co[26] = ci[6];
        co[24] = ci[7];
        co[1] = ci[8];
        co[5] = ci[9];
        co[7] = ci[10];
        co[3] = ci[11];
        co[19] = ci[12];
        co[23] = ci[13];
        co[25] = ci[14];
        co[21] = ci[15];
        co[9] = ci[16];
        co[11] = ci[17];
        co[17] = ci[18];
        co[15] = ci[19];
        co[12] = ci[20];
        co[14] = ci[21];
        co[10] = ci[22];
        co[16] = ci[23];
        co[4] = ci[24];
        co[22] = ci[25];
        co[13] = ci[26];
    }
    return co;
}

std::string pntMesh::getBlockName(int id) const
{
    if (id>=numBlocks || id<0) 
    {
      std::cerr << "Invalid Block Id \n";
      throw;
    }
    return elmBlks[id].regionName;
}

elementType pntMesh::getBlockElmType(int id) const
{
    if (id>=numBlocks || id<0) 
    {
      std::cerr << "Invalid Block Id \n";
      throw;
    }
    return elmBlks[id].eTpe;
}

elementType pntMesh::getElmType(int id) const
{
    if (id>= numElements || id<0)
    {
        std::cerr << "Invalid Element ID \n";
        throw;
    }
    return elmTyp[id];
}

int pntMesh::getElmOrder(int id) const
{
    if (id>= numElements || id<0)
    {
        std::cerr << "Invalid Element ID \n";
        throw;
    }
    return elmOrd[id];
}


VTKCellType pntMesh::getVtkCellTag(elementType et, int order) const
{
    if (et == elementType::TRIANGLE)
      switch(order)
      {
        case 1: return VTK_TRIANGLE;
        case 2: return VTK_QUADRATIC_TRIANGLE;
        default: return VTK_HIGHER_ORDER_TRIANGLE;
      }

    if (et == elementType::QUADRILATERAL)
      switch(order)
      {
        case 1: return VTK_QUAD;
        case 2: return VTK_QUADRATIC_QUAD;
        default: return VTK_HIGHER_ORDER_QUAD;
      }


    if (et == elementType::TETRAHEDRON)
      switch(order)
      {
        case 1: return VTK_TETRA;
        case 2: return VTK_QUADRATIC_TETRA;
        default: return VTK_HIGHER_ORDER_TETRAHEDRON;
      }

    if (et == elementType::HEXAGON || et == elementType::BRICK 
            || et == elementType::LAGRANGE_BRICK)
      switch(order)
      {
        case 1: return VTK_HEXAHEDRON;
        case 2: return VTK_QUADRATIC_HEXAHEDRON;
        default: return VTK_HIGHER_ORDER_HEXAHEDRON;
      }

    if (et == elementType::PRISMATIC)
      switch(order)
      {
        case 1: return VTK_WEDGE;
        case 2: return VTK_QUADRATIC_WEDGE;
        default: return VTK_HIGHER_ORDER_WEDGE;
      }

    std::cerr << "Unknown element type " 
              << et 
              << " order "
              << order
              << std::endl;
    throw;
}

// writes PNT mesh data into file
void pntMesh::write(std::string fname) const
{
  std::ofstream outputStream(fname.c_str());
  if(!outputStream.good()) 
  {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }

  // ---------  writing pntmesh header ----------- //
  outputStream << std::setw(10) << numVertices 
               << std::setw(10) << numElements
               << std::setw(10) << numDimensions
               << std::setw(10) << numBlocks
               << std::setw(10) << numSurfaces
               << std::setw(10) << numSurfInternal
               << std::setw(10) << numSurfBoundary
               << std::setw(10) << 0
               << std::setw(10) << 0
               << std::setw(10) << 0
               << std::endl;

  // ---------  writing node coords ----------- //
  int lb = 0;
  for (int id=0; id<numDimensions; id++)
  {
    for (int iv=0; iv<numVertices; iv++)
    {
      outputStream << std::setw(15)
                   << std::scientific 
                   << std::setprecision(8)
                   << pntCrds[iv][id];
      if (++lb == 8) 
      {
        outputStream << std::endl;
        lb = 0;
      }
    }
  }
  outputStream << std::endl;
  
  // ---------  writing element blocks ----------- //
  for (int ib=0; ib<numBlocks; ib++)
  {
    // block header
    outputStream << std::setw(10)
                 << elmBlks[ib].numElementsInBlock
                 << std::setw(10)
                 << elmBlks[ib].numBoundarySurfacesInBlock
                 << std::setw(10)
                 << elmBlks[ib].nodesPerElement
                 << std::setw(10)
                 << 0
                 << std::setw(10)
                 << elmBlks[ib].ordIntrp
                 << std::setw(10)
                 << elmBlks[ib].ordEquat
                 << std::endl;
    
    // element type and region name
    outputStream << std::setw(16)
                 << std::left
                 << elmTypeStr(elmBlks[ib].eTpe)
                 << std::setw(16)
                 << std::left
                 << elmBlks[ib].regionName
                 << std::endl;

    // element connectivity
    lb = 1;
    for (int ie=0; ie<elmBlks[ib].numElementsInBlock; ie++)
    {
       std::vector<int> econn;
       econn = elmBlks[ib].eConn[ie];
       econn = getPntConn(econn, elmBlks[ib].eTpe, elmBlks[ib].ordIntrp);
       for (int im=0; im<econn.size(); im++, lb++)
       {
         outputStream << std::setw(10)
                      << std::right
                      << econn[im]+1;
         if (lb == 12)
         {
           lb = 0;
           outputStream << std::endl;
         }
       }
    }
    if (lb != 1) outputStream << std::endl;

    // BC tags
    lb = 1;
    for (auto it=elmBlks[ib].srfBCTag.begin(); 
              it!=elmBlks[ib].srfBCTag.end(); it++, lb++)
    {
      outputStream << std::setw(16)
                   << std::left
                   << bcTagStr(*it);
      if (lb == 7)
      {
        lb = 0;
        outputStream << std::endl;
      }
    }
    if (lb != 1) outputStream << std::endl;

    // element number/ref number pairs
    lb = 1;
    for (auto it=elmBlks[ib].srfBCEleRef.begin();
              it!=elmBlks[ib].srfBCEleRef.end();
              it++, lb++)
    {
      outputStream << std::setw(10)
                   << std::right
                   << *it;
      if (lb == 12)
      {
        lb = 0;
        outputStream << std::endl;
      }
    }
    if (lb != 1) outputStream << std::endl;

    // adjacancy header
    outputStream << std::setw(10)
                 << elmBlks[ib].numElementsInBlock
                 << std::setw(10)
                 << elmBlks[ib].numSurfPerEleInBlock
                 << std::endl;

    // global surface array
    lb = 1;
    for (auto it=elmBlks[ib].glbSrfId.begin();
              it!=elmBlks[ib].glbSrfId.end();
              it++,lb++)
    {
      outputStream << std::setw(10)
                   << std::right
                   << *it;
      if (lb == 12)
      {
        lb = 0;
        outputStream << std::endl;
      }
    }
    if (lb != 1) outputStream << std::endl;

    // adjacancy block array
    lb = 1;
    for (auto it=elmBlks[ib].adjBlkId.begin();
              it!=elmBlks[ib].adjBlkId.end();
              it++,lb++)
    {
      outputStream << std::setw(10)
                   << std::right
                   << *it;
      if (lb == 12)
      {
        lb = 0;
        outputStream << std::endl;
      }
    }
    if (lb != 1) outputStream << std::endl;

    // adjacancy element array
    lb = 1;
    for (auto it=elmBlks[ib].adjElmId.begin();
              it!=elmBlks[ib].adjElmId.end();
              it++,lb++)
    {
      outputStream << std::setw(10)
                   << std::right
                   << *it;
      if (lb == 12)
      {
        lb = 0;
        outputStream << std::endl;
      }
    }
    if (lb != 1) outputStream << std::endl;

    // adjacancy reference array
    lb = 1;
    for (auto it=elmBlks[ib].adjRefId.begin();
              it!=elmBlks[ib].adjRefId.end();
              it++,lb++)
    {
      outputStream << std::setw(10)
                   << std::right
                   << *it;
      if (lb == 12)
      {
        lb = 0;
        outputStream << std::endl;
      }
    }
    if (lb != 1) outputStream << std::endl;
  }
  
  // closing the stream
  outputStream.close();

}


void pntMesh::pntPopulate(const meshBase* imb) 
{
  // acquiring dataset
  vtkSmartPointer<vtkDataSet> ds = imb->getDataSet();
  if (!ds)
  {
    std::cerr << "No dataset is associated to the meshbase.\n";
    exit(1);
  }
  
  // loop through cells and obtain different quanitities 
  // needed
  int srfId = 0;
  int nCl = ds->GetNumberOfCells();
  elmSrfId.resize(nCl);
  for (int ic=0; ic<nCl; ic++)
  {
    vtkCell* vc = ds->GetCell(ic); 
    if (vc->GetCellDimension() == 2)
    {
      //for 2D cells
      int ne = vc->GetNumberOfEdges();
      for (int ie=0; ie<ne; ie++)
      {
        numSurfInternal++;
        vtkCell* ve = vc->GetEdge(ie);
        vtkIdList* pidl = ve->GetPointIds();
        
        // edge connectivity
        std::set<int> econn;
        std::map<std::set<int>,int >::const_iterator itSrfId;
        std::pair<int,int> adjPair;
        std::vector<int> adjCellId;
        for (int ipt=0; ipt<pidl->GetNumberOfIds(); ipt++)
          econn.insert(pidl->GetId(ipt));
        auto ret1 = connSet.insert(econn);
        bool isNew = ret1.second;
        if (isNew) 
        {
          // global surface connectivity
          surfConnToId[econn] = srfId;
          surfIdToConn[srfId] = econn;
          // global reference number
          adjPair.first = ic;
          adjPair.second = ie+1;
        } 
        else 
        {
          // if edge is already accounted update reference
          // number for its second neighbor cell
          itSrfId = surfConnToId.find(econn);
          if (itSrfId != surfConnToId.end())
          {
            adjPair.first = ic;
            adjPair.second = ie+1;
          } 
          else 
          {
            std::cerr << "Found a repeated surface connectivity without id.";
            exit(1);
          }
        }

        // getting neighbour list
        bool isBndrSrf = false;
        adjCellId.push_back(ic+1);
        vtkSmartPointer<vtkIdList> cidl = vtkSmartPointer<vtkIdList>::New();
        ds->GetCellNeighbors(ic, pidl, cidl);
        if (cidl->GetNumberOfIds() == 0)
        {
          isBndrSrf = true;
          adjPair.second = 0;
          adjCellId.push_back(0);
          numSurfBoundary++;
        }

        // updating adjacency information
        if (isNew) {
          // element global surface id
          elmSrfId[ic].push_back(srfId);
          surfOnBndr.push_back(isBndrSrf);
          surfAdjRefNum[srfId].push_back(adjPair);
          surfAdjElmNum[srfId].insert(surfAdjElmNum[srfId].end(),  
                                      adjCellId.begin(), 
                                      adjCellId.end() );
        } 
        else 
        {
          // element global surface id
          elmSrfId[ic].push_back(itSrfId->second);
          surfAdjRefNum[itSrfId->second].push_back(adjPair);
          surfAdjElmNum[itSrfId->second].insert(surfAdjElmNum[itSrfId->second].end(),  
                                      adjCellId.begin(), 
                                      adjCellId.end() );
        }

        // finally increment surface id number if needed
        if (isNew) srfId++;
      }
    } 
    else if (vc->GetCellDimension() == 3)
    {
      //for 3D cells
      int nfc = vc->GetNumberOfFaces();
      for (int ifc=0; ifc<nfc; ifc++)
      {
        numSurfInternal++;
        vtkCell* vf = vc->GetFace(ifc);
        vtkIdList* pidl = vf->GetPointIds();

        // surface connectivity
        std::set<int> econn;
        std::map<std::set<int>,int >::const_iterator itSrfId;
        std::pair<int,int> adjPair;
        std::vector<int> adjCellId;
        for (int ipt=0; ipt<pidl->GetNumberOfIds(); ipt++)
          econn.insert(pidl->GetId(ipt));
        auto ret1 = connSet.insert(econn);
        bool isNew = ret1.second;
        if (isNew) 
        {
          // global surface connectivity
          surfConnToId[econn] = srfId;
          surfIdToConn[srfId] = econn;
          // global reference number
          adjPair.first = ic;
          adjPair.second = ifc+1;
        } 
        else 
        {
          // if edge is already accounted update reference
          // number for its second neighbor cell
          itSrfId = surfConnToId.find(econn);
          if (itSrfId != surfConnToId.end())
          {
            adjPair.first = ic;
            adjPair.second = ifc+1;
          } 
          else 
          {
            std::cerr << "Found a repeated surface connectivity without id.";
            exit(1);
          }
        }



        // getting neighbour list
        bool isBndrSrf = false;
        adjCellId.push_back(ic+1);
        vtkSmartPointer<vtkIdList> cidl = vtkSmartPointer<vtkIdList>::New();
        ds->GetCellNeighbors(ic, pidl, cidl);
        if (cidl->GetNumberOfIds() == 0)
        {
          isBndrSrf = true;
          adjPair.second = 0;
          adjCellId.push_back(0);
          numSurfBoundary++;
        }
        
        // updating adjacency information
        if (isNew) {
          // element global surface id
          elmSrfId[ic].push_back(srfId);
          surfOnBndr.push_back(isBndrSrf);
          surfAdjRefNum[srfId].push_back(adjPair);
          surfAdjElmNum[srfId].insert(surfAdjElmNum[srfId].end(),  
                                      adjCellId.begin(), 
                                      adjCellId.end() );
        } 
        else 
        {
          // element global surface id
          elmSrfId[ic].push_back(itSrfId->second);
          surfAdjRefNum[itSrfId->second].push_back(adjPair);
          surfAdjElmNum[itSrfId->second].insert(surfAdjElmNum[itSrfId->second].end(),  
                                      adjCellId.begin(), 
                                      adjCellId.end() );
        }
        // finally increment surface id number if needed
        if (isNew) srfId++;
      }
    }
  }
  numSurfInternal+=numSurfBoundary;
  numSurfInternal/=2;
  numSurfInternal-=numSurfBoundary;
  numSurfaces = numSurfInternal + numSurfBoundary; 
  std::cout << "Number of boundary surfaces (edges) = " << numSurfBoundary << std::endl;
  std::cout << "Number of internal surfaces (edges) = " << numSurfInternal << std::endl;
  std::cout << "Total number of surfaces (edges) = " << numSurfaces << std::endl;
  std::cout << "Srf ID = " << srfId << std::endl;

  // testing
  /*
  for (auto itt1 = surfAdjElmNum.begin(); itt1 != surfAdjElmNum.end(); itt1++)
  {
    std::cout << "/////////\n";
    for (auto itt2 = (itt1->second).begin(); itt2 != (itt1->second).end(); itt2++)
      std::cout << *itt2 << std::endl;
  }
  */

  /*
  for (auto it = elmSrfId.begin(); it!=elmSrfId.end(); it++)
    std::cout << "Number of element surfaces = " << it->size() << std::endl;
  */

  // preparing element blocks
  for (int iBlk=0; iBlk<numBlocks; iBlk++)
    updElmBlk(iBlk);
  std::cout << "Finished processing blocks." << std::endl; 
}

void pntMesh::updElmBlk(int blkId)
{
  std::cout << "Working on block " << blkId << std::endl;
  elmBlks[blkId].numBoundarySurfacesInBlock = 0;
  // block's first element global index
  int frstElmGlbIdx = elmBlks[blkId].elmIds[0];
  // TODO: Removing the first element to avoid additional tag
  surfaceBCTag sbc = elmBlks[blkId].srfBCTag[0];
  elmBlks[blkId].srfBCTag.pop_back();
  // looping through elements in block
  for (auto ie=elmBlks[blkId].elmIds.begin();
            ie!=elmBlks[blkId].elmIds.end();
            ie++)
  {
    // connectivity
    std::cout << "Working on element " << *ie << std::endl;
    elmBlks[blkId].eConn.push_back(elmConn[*ie]);
    
    // surface related calculations
    // looping through element surfaces
    int srfRefId = 0;
    for (auto is=elmSrfId[*ie].begin();
              is!=elmSrfId[*ie].end();
              is++)
    {
      // incrementing reference id
      srfRefId++;
      // adjacent element information
      int adjElmId = surfAdjRefNum[*is][0].first == *ie ? 
          surfAdjRefNum[*is][1].first : surfAdjRefNum[*is][0].first;
      int adjRefId = surfAdjRefNum[*is][0].first == *ie ? 
          surfAdjRefNum[*is][1].second : surfAdjRefNum[*is][0].second;

      if (surfOnBndr[*is])
      {
        // number of boundary surfaces
        elmBlks[blkId].numBoundarySurfacesInBlock++;
        // TODO: duplicating surface boundary condition tag 
        // with index 0 for now
        elmBlks[blkId].srfBCTag.push_back(sbc);
        // element's boundary surface reference number
        elmBlks[blkId].srfBCEleRef.push_back(*ie+1-frstElmGlbIdx); // local element indx
        elmBlks[blkId].srfBCEleRef.push_back(srfRefId); // ref surface number
        // adjacancy information
        elmBlks[blkId].adjElmId.push_back(0);
        elmBlks[blkId].adjBlkId.push_back(0);
        elmBlks[blkId].adjRefId.push_back(0);
      } 
      else
      {
        // adjacancy information
        elmBlks[blkId].adjElmId.push_back(elmLocalId[adjElmId]+1);
        elmBlks[blkId].adjBlkId.push_back(elmBlkId[adjElmId]+1);
        elmBlks[blkId].adjRefId.push_back(adjRefId);
      }

      // global surface id
      elmBlks[blkId].glbSrfId.push_back(*is+1);

    }

  }
}
