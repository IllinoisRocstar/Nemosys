// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_DISTANCEFUNCTION
#define _H_DISTANCEFUNCTION

#include "MeshDataBaseInterface.h"
#include "DistanceToPoints.h"
#include <set>
#include <map>

namespace MAd {

  // -------------------------------------------------------------------
  // This class computes and stores:
  //   - the distance to a set of vertices, edges or faces 
  //     (as attached double),
  //   - the gradient of the distance (as attached pointer),
  //   - the Laplacian of the distance (as attached double).
  //
  // There are two ways to compute the distance: 
  //   - to a set of vertices (inacurate to represent a wall)
  //   - to a set of edges(2D)/faces(3D) (more accurate, more expensive)
  // The choice is governed by the variable 'distToEntities'
  // -------------------------------------------------------------------
  class distanceFunction {

  public:

    distanceFunction(pMesh m, bool _distToEntities);
    ~distanceFunction();

    void computeTree(const std::set<pVertex>&, const std::set<pEntity>&);

    void clearVertexData(pVertex pv) const;

    // Distance

    void   computeAllDistances() const;
    void   computeAllDistancesEDP(const std::set<pGEntity> fixed) const;
    double getDistance(const pVertex pv) const;
    void   clearDistance(const pVertex pv) const;
    void   clearDistances() const;

    double computeDistance(const double xyz[3]) const;
    double computeDistSq  (const double xyz[3]) const;

    void outputDistance(const char * fn) const;

    // Gradient of the distance

    void computeAllDistAndGrad() const; // the most accurate for gradients

    void computeGradientAtVertices(); // not accurate
    void clearGradientAtVertices();
    bool getGradient(const pVertex pv, double grad[3]) const;
    bool getGradientOnEntity(const pEntity entity, 
                             const double xyz[3],
                             double grad[3]) const;
    void attachGradient(pVertex pv, double grad[3]) const;

    void outputGradAtVertices(const char * fn) const;

    // Curvature ( Laplacian of the distance )

    void computeGradientAndCurvature(const std::set<pRegion>& regs); // not accurate...
    void computeGradientAndCurvature2D(const std::set<pFace>& faces);

    void computeCurvature(const std::set<pRegion>& regs);
    void limitCurvature(double maxCurv) const;
    void smoothCurvature(double maxGrad) const;
    void smoothCurvatureDummy(int nbSmoothings=1) const;
    bool getCurvature(const pVertex pv, double *c) const;
    bool getCurvatureOnEntity(const pEntity entity, 
                              const double xyz[3],
                              double *curv) const;
    void clearCurvature() const;
    void attachCurvature(pVertex pv, double curv) const;

    void outputCurvature(const char * fn) const;

  private:

    double computeDistance(const pVertex pv) const;
    void computeGradientInElements() const;
    void clearGradientInElements() const;

    bool getGradientOnEdge(const pEdge edge, 
                           const double xyz[3],
                           double grad[3]) const;
    bool getGradientOnEdgeParam(const pEdge pE, 
                                const double u,
                                double grad[3]) const;
    bool getGradientOnFace(const pFace face, 
                           const double xyz[3],
                           double grad[3]) const;
    bool getGradientOnFaceParam(const pFace face, 
                                const double u[2],
                                double grad[3]) const;
    bool getGradientOnRegion(const pRegion region, 
                             const double xyz[3],
                             double grad[3]) const;
    bool getGradientOnRegionParam(const pRegion region, 
                                  const double u[3],
                                  double grad[3]) const;

  private:

    pMesh mesh;
    pMeshDataId distId, vGradId, rGradId, vCurvId;

    // true:  distance computed to a set of edges (2D) or faces (3D)
    // false: distance computed to a set of vertices
    bool distToEntities;

    int nVert, nEnt;
    int nVE; // number of vertices in a wall entity
    double * xyzV; // coordinates of the vertices sorted by local ids with contiguous xyz
    int * entToV;  // entities vertices (by local id)
    std::multimap<int,int> vToEnt; // vertices entities
    SearchTool * kdSearch;
    std::map<pVertex,int> pvToSearchId;
    
  };

  // -------------------------------------------------------------------
}

#endif
