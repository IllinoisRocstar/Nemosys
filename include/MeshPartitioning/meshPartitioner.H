#ifndef NEMOSYS_MESHPARTITIONER_H_
#define NEMOSYS_MESHPARTITIONER_H_

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>

// other
#ifdef HAVE_METIS
  #include <metis.h>
#endif
#ifdef HAVE_GMSH
#include <MAdLib.h>
#endif

// forward declare cgnsAnalyzer to prevent gmsh name conflicts w/ simmetrix
class cgnsAnalyzer;

enum class MeshType_t {
  MESH_TRI_3,
  MESH_QUAD_4,
  MESH_TETRA_4,
  MESH_HEX_8,
  MESH_MIXED
};

/* Defines class(es) for mesh partitioning */
class NEMOSYS_EXPORT meshPartition
{
  public:
    meshPartition(int pidx,
                  const std::vector<int> &glbElmPartedIdx,
                  const std::vector<int> &glbElmConn,
                  MeshType_t inMshType);

    ~meshPartition() = default;

    std::vector<double> getCrds(const std::vector<double> &crds) const;
    std::vector<int> getConns() const { return partElmConn; }
    std::vector<double> getElmSlns(const std::vector<double> &slns) const;
    std::vector<double> getElmSlnsVec(const std::vector<double> &slns, int nComp) const;

    std::map<int, int> getPartToGlobNodeMap() const { return ndeIdxPartToGlob; }
    std::map<int, int> getPartToGlobElmMap() const { return elmIdxPartToGlob; }

  public:
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


class NEMOSYS_EXPORT meshPartitioner
{
  public:
    // constructors
    // from scratch
    meshPartitioner(int nNde, int nElm,
                    const std::vector<int> &elemConn,
                    MeshType_t meshType) :
        nNde(nNde), nElm(nElm), nPart(0), meshType(meshType)
    {
      elmConn.insert(elmConn.begin(), elemConn.begin(), elemConn.end());
    }

    /* disable -- AEG
    // from CGNS object
    explicit meshPartitioner(cgnsAnalyzer *inCg);
    */
#ifdef HAVE_GMSH
    // from MAdMesh object
    explicit meshPartitioner(MAd::pMesh inMesh);
#endif
    // from meshBase object
    explicit meshPartitioner(const meshBase *inMB);

    // destructor
    ~meshPartitioner()
    {
      if (!meshParts.empty())
        for (int iPart = 0; iPart < nPart; iPart++)
          delete meshParts[iPart];
    };

    // mesh information
    void setNPartition(int nPartition) { nPart = nPartition; }
    int getNPartition() const { return nPart; }
    int partition(int nPartition);
    int partition();

    std::vector<double> getPartedNde() const;
    std::vector<double> getPartedElm() const;

    void setPartedElm(const std::vector<double> &prtElm);

    std::vector<double> getCrds(int iPart,
                                const std::vector<double> &crds) const
    { return meshParts[iPart]->getCrds(crds); }
    std::vector<int> getConns(int iPart) const
    { return meshParts[iPart]->getConns(); }
    int getNNdePart(int iPart) const { return meshParts[iPart]->nNde; }
    int getNElmPart(int iPart) const { return meshParts[iPart]->nElm; }

    std::map<int, int> getPartToGlobNodeMap(int iPart) const;
    std::map<int, int> getPartToGlobElmMap(int iPart) const;
    // solution data

    // for scalar nodal values the same operator acting on
    // coordinates can be used to find solution values for
    // the partition
    std::vector<double> getNdeSlnScalar(int iPart,
                                        const std::vector<double> &slns) const
    { return getCrds(iPart, slns); }

    std::vector<double> getElmSlnScalar(int iPart,
                                        const std::vector<double> &slns) const
    { return meshParts[iPart]->getElmSlns(slns); }
    std::vector<double> getElmSlnVec(int iPart,
                                     const std::vector<double> &slns,
                                     int nComp) const
    { return meshParts[iPart]->getElmSlnsVec(slns, nComp); }

  private:
    void buildPartitions();

  private:
    int nNde;
    int nElm;
    int nPart;
    std::vector<int> elmConnVec;
    std::vector<int> elmConn;
    MeshType_t meshType;
#ifdef HAVE_METIS
    // metis datastructures
    idx_t options[METIS_NOPTIONS];
#endif
    std::vector<int> epart;
    std::vector<int> npart;
    // partition data
    std::vector<meshPartition *> meshParts;
};

#endif  // NEMOSYS_MESHPARTITIONER_H_
