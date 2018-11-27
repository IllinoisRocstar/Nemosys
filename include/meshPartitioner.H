/* Defines class(es) for mesh partitioning */

#ifndef __MESHPARTITIONER_H
#define __MESHPARTITIONER_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>

// other
#include "metis.h"
#include "MAdLib.h"

// forward declare cgnsAnalyzer to prevent gmsh name conflicts w/ simmetrix
class cgnsAnalyzer;
class meshBase;

enum MeshType_t {MESH_TRI_3, MESH_QUAD_4, MESH_TETRA_4, MESH_HEX_8, MESH_MIXED};

class meshPartition{
public:

  meshPartition(int pidx, std::vector<int> glbElmPartedIdx, std::vector<int> glbElmConn, MeshType_t inMshType);
 ~meshPartition(){};

  std::vector<double> getCrds(std::vector<double>& crds);
  std::vector<int> getConns();
  std::vector<double> getElmSlns(std::vector<double>& slns);
  std::vector<double> getElmSlnsVec(std::vector<double>& slns, int nComp);
  
  std::map<int,int> getPartToGlobNodeMap();
  std::map<int,int> getPartToGlobElmMap();
  int nNde;
  int nElm;

private:
  int pIdx;
  int nNdeElm;
  std::vector<int> globNdeIdx;
  std::vector<int> globElmIdx;
  std::vector<int> partElmConn;
  std::map<int, int> ndeIdxGlobToPart;
  std::map<int, int> ndeIdxPartToGlob;
  std::map<int, int> elmIdxGlobToPart;
  std::map<int, int> elmIdxPartToGlob;
  MeshType_t mshType;
};

class meshPartitioner{

public:
  // constructors
  // from scratch
  meshPartitioner(int nNde, int nElm, std::vector<int>& elemConn, MeshType_t meshType) :
  nNde(nNde), nElm(nElm), nPart(0), meshType(meshType)
  {
    elmConn.insert(elmConn.begin(), elemConn.begin(), elemConn.end());     
  };
  // from CGNS object
  meshPartitioner(cgnsAnalyzer* inCg);
  // from MAdMesh object
  meshPartitioner(MAd::pMesh inMesh);
  //// from meshBase object
  meshPartitioner(const meshBase* inMB);
  // destructor
  ~meshPartitioner()
  {
    if (!meshParts.empty())
      for (int iPart = 0; iPart<nPart; iPart++)
        delete meshParts[iPart];
  };

  // mesh information
  void setNPartition(int nPartition);
  int getNPartition();
  int partition(int nPartition);
  int partition();
  std::vector<double> getPartedNde();
  std::vector<double> getPartedElm();
  void setPartedElm(std::vector<double>& prtElm);
  std::vector<double> getCrds(int iPart, std::vector<double> crds);
  std::vector<int> getConns(int iPart);
  int getNNdePart(int iPart);
  int getNElmPart(int iPart);
  
  std::map<int,int> getPartToGlobNodeMap(int iPart);
  std::map<int,int> getPartToGlobElmMap(int iPart);
  // solution data
  std::vector<double> getNdeSlnScalar(int iPart, std::vector<double>& slns);
  std::vector<double> getElmSlnScalar(int iPart, std::vector<double>& slns);
  std::vector<double> getElmSlnVec(int iPart, std::vector<double>& slns, int nComp);

private:
  void buildPartitions();
   
private:
  int nNde;
  int nElm;
  int nPart;
  std::vector<int> elmConnVec;
  std::vector<int> elmConn;
  MeshType_t meshType;
  // metis datastructures
  idx_t options[METIS_NOPTIONS];
  std::vector<int> epart;
  std::vector<int> npart;
  // partition data
  std::vector<meshPartition* > meshParts;
};
#endif
