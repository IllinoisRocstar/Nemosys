// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef DISCRETE_FACE_H
#define DISCRETE_FACE_H

#include <algorithm>
#include "gmsh/GmshConfig.h"
#include "gmsh/GModel.h"
#include "gmsh/GFace.h"
#include "gmsh/discreteVertex.h"
#include "gmsh/discreteEdge.h"
#include "gmsh/MTriangle.h"
#include "gmsh/MElementCut.h"
#include "gmsh/MEdge.h"
#include "gmsh/MLine.h"
#include "gmsh/rtree.h"

#if defined(MAD_HAVE_HXT)

extern "C" {
#include "hxt_mesh.h"
#include "hxt_parametrization.h"
#include "hxt_linear_system.h"
#include "hxt_curvature.h"
}

class MElementOctree;

class hxt_reparam_surf {
public:
  MElementOctree *oct;
  mutable RTree<std::pair<MTriangle *, MTriangle *> *, double, 3> rtree3d;
  std::vector<MVertex> v2d;
  std::vector<MVertex> v3d;
  std::vector<MTriangle> t2d;
  std::vector<MTriangle> t3d;
  std::vector<GEdge *> bnd;
  std::vector<GEdge *> emb;
  hxt_reparam_surf() : oct(NULL) {}
  ~hxt_reparam_surf();
};

#endif

class discreteFace : public GFace {
private:
  bool _checkAndFixOrientation();
#if defined(MAD_HAVE_HXT)
  int _currentParametrization;
  std::vector<hxt_reparam_surf> _parametrizations;
  HXTStatus _reparametrizeThroughHxt();
  bool _computeTopologyOfPartition(int nbColors, int *colors, int *nNodes,
                                   int *nodes, double *uv,
                                   std::vector<MVertex *> &c2v,
                                   std::vector<std::vector<MEdge> > &boundaries);
#endif
public:
  discreteFace(GModel *model, int num);
  virtual ~discreteFace() {}
  using GFace::point;
  GPoint point(double par1, double par2) const;
  SPoint2 parFromPoint(const SPoint3 &p, bool onSurface = true) const;
  GPoint closestPoint(const SPoint3 &queryPoint, double maxDistance,
                      SVector3 *normal = NULL) const;
  GPoint closestPoint(const SPoint3 &queryPoint,
                      const double initialGuess[2]) const;
  SVector3 normal(const SPoint2 &param) const;
  double curvatureMax(const SPoint2 &param) const;
  double curvatures(const SPoint2 &param, SVector3 &dirMax, SVector3 &dirMin,
                    double &curvMax, double &curvMin) const;
  GEntity::GeomType geomType() const { return DiscreteSurface; }
  virtual Pair<SVector3, SVector3> firstDer(const SPoint2 &param) const;
  virtual void secondDer(const SPoint2 &param, SVector3 &dudu, SVector3 &dvdv,
                         SVector3 &dudv) const;
  void createGeometry();
  virtual bool haveParametrization()
  {
#if defined(MAD_HAVE_HXT)
    return !_parametrizations.empty();
#else
    return false;
#endif
  }
  virtual void mesh(bool verbose);
  void setBoundEdges(const std::vector<int> &tagEdges);
  void setBoundEdges(const std::vector<int> &tagEdges,
                     const std::vector<int> &signEdges);
  int trianglePosition(double par1, double par2, double &u, double &v) const;
  GPoint intersectionWithCircle(const SVector3 &n1, const SVector3 &n2,
                                const SVector3 &p, const double &R,
                                double uv[2]);
};

#endif
