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
// Authors: Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef _H_MODELINTERFACE
#define _H_MODELINTERFACE

#include "madlib_export.h"

#include <string>
#include <vector>

namespace MAd {

  // -------------------------------------------------------------------
#ifdef _HAVE_GMSH_
  typedef class GmshGEntity MAdGEntity;
  typedef class GmshGEntityLessThan MAdGEntityLessThan;
  typedef class GmshGRegion MAdGRegion;
  typedef class GmshGFace   MAdGFace;
  typedef class GmshGEdge   MAdGEdge;
  typedef class GmshGVertex MAdGVertex;
#else
  typedef class NullGEntity MAdGEntity;
  typedef class NullGEntityLessThan MAdGEntityLessThan;
  typedef class NullGRegion MAdGRegion;
  typedef class NullGFace   MAdGFace;
  typedef class NullGEdge   MAdGEdge;
  typedef class NullGVertex MAdGVertex;
#endif
  
  typedef MAdGEntity * pGEntity;
  typedef MAdGRegion * pGRegion;
  typedef MAdGFace   * pGFace;
  typedef MAdGEdge   * pGEdge;
  typedef MAdGVertex * pGVertex;

  typedef class MAdModel * pGModel;

  typedef class GM_RegionIterator * GRIter;
  typedef class GM_FaceIterator   * GFIter;
  typedef class GM_EdgeIterator   * GEIter;
  typedef class GM_VertexIterator * GVIter;

  typedef class PGList * pPGList;

  // -------------------------------------------------------------------

  // --- Model operators ---

  // create an empty model
  void MADLIB_EXPORT GM_create(pGModel * model, std::string name="");

  // delete the model
  void MADLIB_EXPORT GM_delete(pGModel model);

  // read a model, guess file format from extension
  int  MADLIB_EXPORT GM_read(pGModel model, const std::string name);

  // read particular file formats
  int  MADLIB_EXPORT GM_readFromMSH(pGModel model, const std::string name);
  int  MADLIB_EXPORT GM_readFromGEO(pGModel model, const std::string name);
  int  MADLIB_EXPORT GM_readFromSTEP(pGModel model, const std::string name);
  int  MADLIB_EXPORT GM_readFromBREP(pGModel model, const std::string name);
  int  MADLIB_EXPORT GM_readFromIGES(pGModel model, const std::string name);

  // void GM_diagnostics(const pGModel model);

  // Returns true if physical tags are used in the model
  bool MADLIB_EXPORT GM_physical(const pGModel model);

  // Find the entity with the given tag. 
  // Create it if it doesn't exist.
  pGEntity MADLIB_EXPORT GM_entityByTag(const pGModel model, int type, int tag);
  pGRegion MADLIB_EXPORT GM_regionByTag(const pGModel model, int tag);
  pGFace   MADLIB_EXPORT GM_faceByTag  (const pGModel model, int tag);
  pGEdge   MADLIB_EXPORT GM_edgeByTag  (const pGModel model, int tag);
  pGVertex MADLIB_EXPORT GM_vertexByTag(const pGModel model, int tag);

  int MADLIB_EXPORT GM_numVertices(const pGModel);
  int MADLIB_EXPORT GM_numEdges(const pGModel);
  int MADLIB_EXPORT GM_numFaces(const pGModel);
  int MADLIB_EXPORT GM_numRegions(const pGModel);

  // --- Iterators ---

  GRIter MADLIB_EXPORT GM_regionIter(pGModel);
  GFIter MADLIB_EXPORT GM_faceIter(pGModel);
  GEIter MADLIB_EXPORT GM_edgeIter(pGModel);
  GVIter MADLIB_EXPORT GM_vertexIter(pGModel);

  pGRegion MADLIB_EXPORT GRIter_next(GRIter);
  pGFace   MADLIB_EXPORT GFIter_next(GFIter);
  pGEdge   MADLIB_EXPORT GEIter_next(GEIter);
  pGVertex MADLIB_EXPORT GVIter_next(GVIter);

  void MADLIB_EXPORT GRIter_delete(GRIter);
  void MADLIB_EXPORT GFIter_delete(GFIter);
  void MADLIB_EXPORT GEIter_delete(GEIter);
  void MADLIB_EXPORT GVIter_delete(GVIter);

  void MADLIB_EXPORT GRIter_reset(GRIter);
  void MADLIB_EXPORT GFIter_reset(GFIter);
  void MADLIB_EXPORT GEIter_reset(GEIter);
  void MADLIB_EXPORT GVIter_reset(GVIter);

  // --- Entity operators ---

  int  MADLIB_EXPORT GEN_tag(const pGEntity);
  int  MADLIB_EXPORT GEN_type(const pGEntity);
  void MADLIB_EXPORT GEN_setPhysical(pGEntity, int dim, int tag);
  int  MADLIB_EXPORT GEN_physTag(const pGEntity);
  int  MADLIB_EXPORT GEN_physDim(const pGEntity);

#ifdef _HAVE_GMSH_
  std::vector<pGEntity> MADLIB_EXPORT GEN_closure(const pGEntity);

  // --- Region operators ---

  std::vector<pGFace> MADLIB_EXPORT GR_faces(const pGRegion);

  // --- Face operators ---

  int MADLIB_EXPORT GF_numRegions(const pGFace);
  std::vector<pGEdge> MADLIB_EXPORT GF_edges(const pGFace);
  bool MADLIB_EXPORT GF_getParams(const pGFace, const double[3], double[2]);
  void MADLIB_EXPORT GF_closestPoint(const pGFace, const double[3],
                                     const double[2], double[3]);
  void MADLIB_EXPORT GF_xyz(const pGFace, double, double, double[3]);
  double MADLIB_EXPORT GF_curvatureDiv(const pGFace, const double[2],
                                       double cMaxBound);
  double MADLIB_EXPORT GF_curvatures(const pGFace, const double[2],
                                     double dirMax[3], double dirMin[3],
                                     double *curvMax, double *curvMin, 
                                     double cMaxBound);
  void MADLIB_EXPORT GF_centerOnGeodesic(const pGFace face, double t,
                                         const double e[2][2], double c[2]);

  // --- Edge operators ---

  std::vector<pGVertex> MADLIB_EXPORT GE_vertices(const pGEdge);
  void MADLIB_EXPORT GE_closestPoint(const pGEdge, const double[3], double[3]);
  void MADLIB_EXPORT GE_xyz(const pGEdge, double, double[3]);
  void MADLIB_EXPORT GE_reparamOnFace(const pGEdge, const pGFace,
                                      double, double[2], double uClose[2]=NULL);
  bool MADLIB_EXPORT GE_isSeam(const pGEdge, const pGFace);
  double MADLIB_EXPORT GE_curvature(const pGEdge, double, double);

  // --- Vertex operators ---

  std::vector<pGEdge> MADLIB_EXPORT GV_edges(const pGVertex);

  void MADLIB_EXPORT GV_reparamOnFace(const pGVertex, const pGFace, double [2],
                                      double uClose[2]=NULL);
  void MADLIB_EXPORT GV_reparamOnEdge(const pGVertex, const pGEdge, double *,
                                      double uClose=-1.);
  bool MADLIB_EXPORT GV_isOnSeam(const pGVertex, const pGFace);
#else
  void MADLIB_EXPORT GF_centerOnGeodesic(const pGFace face, double t,
                                         const double e[2][2], double c[2]);
#endif

  // --- pGList operators ---

  pPGList  MADLIB_EXPORT PGList_new();
  pPGList  MADLIB_EXPORT PGList_allocate();
  void     MADLIB_EXPORT PGList_delete          (pPGList);
  void     MADLIB_EXPORT PGList_deallocate      (pPGList);
  void     MADLIB_EXPORT PGList_clear           (pPGList);
  pPGList  MADLIB_EXPORT PGList_appPGListUnique (pPGList, pPGList source);
  pPGList  MADLIB_EXPORT PGList_appUnique       (pPGList, pGEntity);
  pPGList  MADLIB_EXPORT PGList_append          (pPGList, pGEntity);
  int      MADLIB_EXPORT PGList_size            (pPGList); 
  pGEntity MADLIB_EXPORT PGList_item            (pPGList, int n);
  pGEntity MADLIB_EXPORT PGList_next            (pPGList, void ** restart);
  int      MADLIB_EXPORT PGList_inList          (pPGList, pGEntity);
  void     MADLIB_EXPORT PGList_remItem         (pPGList, pGEntity);

  // -------------------------------------------------------------------

}

#endif
