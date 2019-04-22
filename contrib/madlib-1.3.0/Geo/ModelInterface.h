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
  void GM_create(pGModel * model, std::string name="");

  // delete the model
  void GM_delete(pGModel model);

  // read a model, guess file format from extension
  int  GM_read(pGModel model, const std::string name);

  // read particular file formats
  int  GM_readFromMSH(pGModel model, const std::string name);
  int  GM_readFromGEO(pGModel model, const std::string name);
  int  GM_readFromSTEP(pGModel model, const std::string name);
  int  GM_readFromBREP(pGModel model, const std::string name);
  int  GM_readFromIGES(pGModel model, const std::string name);

  // void GM_diagnostics(const pGModel model);

  // Returns true if physical tags are used in the model
  bool GM_physical(const pGModel model);

  // Find the entity with the given tag. 
  // Create it if it doesn't exist.
  pGEntity GM_entityByTag(const pGModel model, int type, int tag);
  pGRegion GM_regionByTag(const pGModel model, int tag);
  pGFace   GM_faceByTag  (const pGModel model, int tag);
  pGEdge   GM_edgeByTag  (const pGModel model, int tag);
  pGVertex GM_vertexByTag(const pGModel model, int tag);

  int GM_numVertices(const pGModel);
  int GM_numEdges(const pGModel);
  int GM_numFaces(const pGModel);
  int GM_numRegions(const pGModel);

  // --- Iterators ---

  GRIter GM_regionIter(pGModel);
  GFIter GM_faceIter(pGModel);
  GEIter GM_edgeIter(pGModel);
  GVIter GM_vertexIter(pGModel);

  pGRegion GRIter_next(GRIter);
  pGFace GFIter_next(GFIter);
  pGEdge GEIter_next(GEIter);
  pGVertex GVIter_next(GVIter);

  void GRIter_delete(GRIter);
  void GFIter_delete(GFIter);
  void GEIter_delete(GEIter);
  void GVIter_delete(GVIter);

  void GRIter_reset(GRIter);
  void GFIter_reset(GFIter);
  void GEIter_reset(GEIter);
  void GVIter_reset(GVIter);

  // --- Entity operators ---

  int GEN_tag(const pGEntity);
  int GEN_type(const pGEntity);
  void GEN_setPhysical(pGEntity, int dim, int tag);
  int GEN_physTag(const pGEntity);
  int GEN_physDim(const pGEntity);

#ifdef _HAVE_GMSH_
  std::vector<pGEntity> GEN_closure(const pGEntity);

  // --- Region operators ---

  std::vector<pGFace> GR_faces(const pGRegion);

  // --- Face operators ---

  int GF_numRegions(const pGFace);
  std::vector<pGEdge> GF_edges(const pGFace);
  bool GF_getParams(const pGFace, const double[3], double[2]);
  void GF_closestPoint(const pGFace, const double[3],
                       const double[2], double[3]);
  void GF_xyz(const pGFace, double, double, double[3]);
  double GF_curvatureDiv(const pGFace, const double[2], 
                         double cMaxBound);
  double GF_curvatures(const pGFace, const double[2],
                       double dirMax[3], double dirMin[3],
                       double *curvMax, double *curvMin, 
                       double cMaxBound);
  void GF_centerOnGeodesic(const pGFace face, double t,
                           const double e[2][2], double c[2]);

  // --- Edge operators ---

  std::vector<pGVertex> GE_vertices(const pGEdge);
  void GE_closestPoint(const pGEdge, const double[3], double[3]);
  void GE_xyz(const pGEdge, double, double[3]);
  void GE_reparamOnFace(const pGEdge, const pGFace, 
                        double, double[2], double uClose[2]=NULL);
  bool GE_isSeam(const pGEdge, const pGFace);
  double GE_curvature(const pGEdge, double, double);

  // --- Vertex operators ---

  std::vector<pGEdge> GV_edges(const pGVertex);

  void GV_reparamOnFace(const pGVertex, const pGFace, double [2], 
                        double uClose[2]=NULL);
  void GV_reparamOnEdge(const pGVertex, const pGEdge, double *, 
                        double uClose=-1.);
  bool GV_isOnSeam(const pGVertex, const pGFace);
#else
  void GF_centerOnGeodesic(const pGFace face, double t,
                           const double e[2][2], double c[2]);
#endif

  // --- pGList operators ---

  pPGList  PGList_new();
  pPGList  PGList_allocate();
  void     PGList_delete          (pPGList);
  void     PGList_deallocate      (pPGList);
  void     PGList_clear           (pPGList);
  pPGList  PGList_appPGListUnique (pPGList, pPGList source);
  pPGList  PGList_appUnique       (pPGList, pGEntity);
  pPGList  PGList_append          (pPGList, pGEntity);
  int      PGList_size            (pPGList); 
  pGEntity PGList_item            (pPGList, int n);
  pGEntity PGList_next            (pPGList, void ** restart);
  int      PGList_inList          (pPGList, pGEntity);
  void     PGList_remItem         (pPGList, pGEntity);

  // -------------------------------------------------------------------

}

#endif
