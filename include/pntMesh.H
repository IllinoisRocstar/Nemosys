#ifndef pntMesh_H
#define pntMesh_H

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <vtkCellTypes.h>

class meshBase;

namespace PNTMesh
{

typedef int idTyp;

enum elementType {
  BAR,
  QUADRILATERAL,
  TRIANGLE,
  HEXAGON,
  SPHERICAL,
  CYLINDRICAL,
  BRICK,
  LAGRANGE_BRICK,
  TETRAHEDRON,
  HEXPRISM,
  PRISMATIC,
  OTHER
};

enum surfaceBCTag {
  REFLECTIVE,
  VOID
};

struct blockType {
  int numElementsInBlock;
  int numBoundarySurfacesInBlock;
  int nodesPerElement;
  int ordIntrp;
  int ordEquat;
  int numSurfPerEleInBlock;
  elementType eTpe;
  std::string regionName;
  std::vector<int> elmIds;
  std::vector<surfaceBCTag> srfBCTag;
  std::vector<int> srfBCEleRef;
  std::vector<int> glbSrfId;
  std::vector<int> adjBlkId;
  std::vector<int> adjElmId;
  std::vector<int> adjRefId;
  std::vector<std::vector<int> > eConn;
};

typedef std::vector<blockType> BlockMap;

VTKCellType p2vEMap (elementType et);
elementType v2pEMap (VTKCellType vt);
surfaceBCTag bcTagNum(std::string& tag);
std::string bcTagStr(int tag);
elementType elmTypeNum (std::string tag);
std::string elmTypeStr (elementType tag);
int elmNumNde (elementType et, int order);
int elmNumSrf (elementType et);

class pntMesh
{

    public:
      pntMesh();
      pntMesh(std::string ifname);
      pntMesh(const meshBase* imb, int dim, int nBlk, BlockMap& elmBlkMap);
      ~pntMesh() {};

    public:
      bool isCompatible() const {return isSupported;}
      int getNumberOfPoints() const;
      int getNumberOfCells() const;
      int getNumberOfBlocks() const {return numBlocks;}

      std::vector<double> getPointCrd(int id) const;
      std::vector<int> getElmConn(int id) const;
      std::vector<int> getElmConn(int id, VTKCellType vct) const;
      std::vector<int> getPntConn(std::vector<int>& ci, elementType et, int eo) const;
      std::string getBlockName(int id) const;
      elementType getBlockElmType(int id) const;

      elementType getElmType(int id) const;
      int getElmOrder(int id) const;
      VTKCellType getVtkCellTag(elementType et, int order) const;

      void write(std::string fname) const;

   private:
      // populate quanitites needed to write
      // pntmesh file
      void pntPopulate(const meshBase* imb);
      void updElmBlk(int blkId);

    private:
      std::string ifname;
      int numVertices; 
      int numElements; 
      int numDimensions;
      int numBlocks; 
      int numSurfaces; 
      int numSurfInternal; 
      int numSurfBoundary;
      std::vector<std::vector<double> > pntCrds;
      std::vector<std::vector<int> > elmConn;    // connvectivity
      std::vector<elementType> elmTyp;           // type
      std::vector<int> elmOrd;                   // order
      std::vector<int> elmBlkId;                 // element block number
      std::vector<int> elmLocalId;               // elements local id indexed by global ids
      std::vector<std::vector<int> > elmSrfId;   // element surface/edge Ids
      std::vector<blockType> elmBlks;

      bool isSupported; // false is non Tri/Tet elements were found

      // pnt topological information
      std::set<std::set<int> > connSet;
      std::map<std::set<int>,int > surfConnToId;
      std::map<int,std::set<int> > surfIdToConn;
      std::vector<bool> surfOnBndr;
      std::map<int,std::vector<std::pair<int,int> > > surfAdjRefNum;
      std::map<int,std::vector<int> > surfAdjElmNum;
};

} // end namespace pntMesh
#endif
