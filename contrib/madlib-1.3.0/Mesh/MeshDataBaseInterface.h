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

#ifndef H_MESHDATABASEINTEFACE
#define H_MESHDATABASEINTEFACE

#include "ModelInterface.h"
#include <iostream>
#include <set>

namespace MAd {

  class MDB_DataExchangerPeriodic;

  typedef unsigned int pMeshDataId;

  typedef class MDB_Mesh  * pMesh;  

  typedef class MDB_MeshEntity * pEntity;
  typedef class MDB_Region     * pRegion;
  typedef class MDB_Face       * pFace;
  typedef class MDB_Edge       * pEdge;
  typedef class MDB_Point      * pVertex;
  typedef class MDB_Point      * pPoint;
  
  typedef class MDB_RegionIter * RIter;
  typedef class MDB_FaceIter   * FIter;
  typedef class MDB_EIter      * EIter;
  typedef class MDB_VIter      * VIter;

  typedef class PList * pPList;

  /********************/
  /*  Mesh Operators  */
  /********************/

  pMesh   M_new(pGModel);
  void    M_delete(pMesh);
  void    M_load(pMesh, const char * filename);
  void    M_clean(pMesh);

  void    M_shrink(pMesh);
  void    M_expand(pMesh);

  void    M_info(const pMesh, std::ostream& out=std::cout);
  pGModel M_model(pMesh);
  int     M_dim(pMesh);

  int     M_edgesMaxOrder(pMesh);
  int     M_facesMaxOrder(pMesh);
  int     M_regionsMaxOrder(pMesh);
  int     M_maxOrder(pMesh);

  bool    M_isParametric(pMesh);

  void    M_writeMsh(const pMesh, const char *name, int version=2, 
                     const int * partitionTable=NULL);
  void    M_writeMshPer(pMesh, const char *name, 
                        MDB_DataExchangerPeriodic &, int version);

#ifdef _HAVE_METIS_
  void    M_Partition(pMesh, int nbParts, const char * filename);
#endif

  int M_numRegions(pMesh);
  int M_numTets(pMesh);
  int M_numHexes(pMesh);
  int M_numPrisms(pMesh);
  int M_numFaces(pMesh);
  int M_numQuads(pMesh);
  int M_numTriangles(pMesh);
  int M_numEdges(pMesh);
  int M_numVertices(pMesh);

  RIter M_regionIter(pMesh);
  FIter M_faceIter(pMesh);
  EIter M_edgeIter(pMesh);
  VIter M_vertexIter(pMesh);

  int M_numClassifiedRegions (pMesh, pGEntity);
  int M_numClassifiedFaces   (pMesh, pGEntity);
  int M_numClassifiedEdges   (pMesh, pGEntity);
  int M_numClassifiedVertices(pMesh, pGEntity);

  RIter   M_classifiedRegionIter(pMesh, pGEntity);
  FIter   M_classifiedFaceIter  (pMesh, pGEntity, int closure);
  EIter   M_classifiedEdgeIter  (pMesh, pGEntity, int closure);
  VIter   M_classifiedVertexIter(pMesh, pGEntity, int closure);
  pVertex M_classifiedVertex    (pMesh, pGVertex);

  void M_classifyEntities(pMesh);

  /* mesh entity creation routines */
  pRegion M_createR(pMesh, int nFace, pFace   faces[], pGEntity pg=NULL);
  pRegion M_createR(pMesh, int nVert, pVertex verts[], pGEntity pg=NULL);
  pRegion M_createR(pMesh, int nVert, int   vertIds[], pGEntity pg=NULL);
  pRegion M_createTet(pMesh, int order, bool serendip, int vertIds[], pGEntity pg=NULL);
  pFace   M_createF(pMesh, int nEdge, pEdge   edges[], pGEntity pg=NULL);
  pFace   M_createF(pMesh, int nVert, pVertex verts[], pGEntity pg=NULL);
  pFace   M_createF(pMesh, int nVert, int   vertIds[], pGEntity pg=NULL);
  pFace   M_createTri(pMesh, int order, int vertIds[], pGEntity pg=NULL);
  pEdge   M_createE(pMesh, pVertex v1, pVertex v2, pGEntity pg=NULL);
  pEdge   M_createE(pMesh, int   vId1, int   vId2, pGEntity pg=NULL);
  pEdge   M_createE(pMesh, pVertex v1, pVertex v2, int order=1, 
                    pVertex pts[]=NULL, pGEntity pg=NULL);
  pVertex M_createV(pMesh, double x, double y, double z, int patch, pGEntity pg=NULL);
  pVertex M_createV2(pMesh, double xyz[3], int patch, pGEntity pg=NULL);
  pVertex M_createVP(pMesh, double x, double y, double z, 
                     double u, double v, int patch, pGEntity pg=NULL);

  /* mesh entity deletion routines */
  void M_removeRegion(pMesh, pRegion);
  void M_removeFace  (pMesh, pFace);
  void M_removeEdge  (pMesh, pEdge);
  void M_removeVertex(pMesh, pVertex);
  void P_delete(pPoint);
  pPoint P_new(void);

  /* extra access to entities */
  pRegion M_region(pMesh, int []);
  pFace   M_face  (pMesh, int []);
  pEdge   M_edge  (pMesh, int, int);
  pVertex M_vertex(pMesh, int);


  /***********************/
  /* Geometric model ops */
  /***********************/
  /* deprecated (rather use ModelInterface) */

  void    M_setGeomFeature    (pMesh, int, pGEntity);
  int     M_numGeomFeatures   (pMesh);
  int     M_geomFeatureId     (pMesh, int ith);
  pPGList M_geomFeature       (pMesh, int id);
  void    M_printGeomFeatures (pMesh);

  /********************/
  /* Entity Iter Ops  */
  /********************/

  pRegion RIter_next  (RIter);
  void    RIter_delete(RIter);
  void    RIter_reset (RIter);
  pFace   FIter_next  (FIter);
  void    FIter_delete(FIter);
  void    FIter_reset (FIter);
  pEdge   EIter_next  (EIter);
  void    EIter_delete(EIter);
  void    EIter_reset (EIter);
  pVertex VIter_next  (VIter);
  void    VIter_delete(VIter);
  void    VIter_reset (VIter);

  /********************/
  /* Entity Operators */
  /********************/
  
  void     EN_info(const pEntity, std::string name="", std::ostream& out=std::cout);
  void     EN_setID(pEntity, int id);
  int      EN_id   (pEntity);
  void     EN_print(pEntity);

  int      EN_type(pEntity); // obsolete: returns dim, use EN_dim() instead
  int      EN_dim(pEntity);
  int      EN_mshTag(pEntity);
  int      EN_whatInType(pEntity);
  pGEntity EN_whatIn(pEntity);
  void     EN_setWhatIn(pEntity,pGEntity);
  
  int      EN_numVertices(const pEntity);
  pPList   EN_vertices(const pEntity);
  int      EN_inClosure(pEntity, pEntity);

  /********************/
  /* Region Operators */
  /********************/

  void     R_info         (const pRegion, std::string name="", std::ostream& out=std::cout);
  void     R_info_quality (const pRegion, std::ostream& out=std::cout);
  void     R_info_topology(const pRegion, std::ostream& out=std::cout);

  pGRegion R_whatIn(pRegion);
  int      R_whatInType(pRegion);
  void     R_setWhatIn(pRegion, pGEntity);

  int      R_numPoints(pRegion);
  pPoint   R_point(pRegion,int);
  int      R_numVertices(pRegion);
  pPList   R_vertices(pRegion);
  pVertex  R_vertex(pRegion, int);
  pVertex  R_fcOpVt(const pRegion, const pFace);
  pPList   R_edges(pRegion);
  pEdge    R_edge(pRegion, int);
  pEdge    R_gtOppEdg(const pRegion, const pEdge);
  int      R_numFaces(pRegion);
  pFace    R_face(pRegion, int);
  pPList   R_faces(pRegion);
  pFace    R_vtOpFc(const pRegion, const pVertex);
  int      R_inClosure(pRegion, pEntity);
  int      R_faceDir(pRegion, int);
  int      R_faceOri(pRegion, int);
  int      R_dirUsingFace(pRegion,pFace);
  int      R_oriUsingFace(pRegion,pFace);

  void     R_coordP1(const pRegion, double [][3]);
  void     R_coord(const pRegion, double [][3]);
  double   R_volume(const pRegion);
  double   R_XYZ_volume (const double [][3]);
  double   R_circumRad(const pRegion);
  double   R_inscrRad(const pRegion);
  bool     R_meanRatioCube(const pRegion, double *);
  bool     R_XYZ_meanRatioCube(const double [][3], double *);
  bool     R_XYZ_isFlat(const double [][3]);
  void     R_linearParams(const pRegion, const double[3], double[3]);
  void     R_center(const pRegion, double [3]);
  void     R_jacobian(const pRegion, double [3][3]);
  double   R_invJacobian(const pRegion, double [3][3]);
  void     R_box(const pRegion, double [3][2]);
  bool     R_inBox(const pRegion, const double [3], double);
  bool     R_contains(const pRegion, const double [3], double);

  /********************/
  /*  Face Operators  */
  /********************/

  void     F_info(const pFace, std::string name="", std::ostream& out=std::cout);

  pFace    F_exist(pEntity, pEntity, pEntity, pEntity);
  pFace    F_exist(pEdge,pVertex);

  pGEntity F_whatIn(pFace);
  int      F_whatInType(pFace);
  void     F_setWhatIn(pFace, pGEntity);

  pPoint   F_point(pFace, int);
  int      F_numPoints(pFace);
  int      F_numPointsTot(pFace);
  int      F_numVertices(pFace);
  pPList   F_vertices(pFace, int dir=1);
  pVertex  F_vertex(pFace, int);
  pEdge    F_vtOpEd(const pFace, const pVertex);
  int      F_numEdges(pFace);
  pPList   F_edges(pFace);
  pEdge    F_edge(pFace, int);
  int      F_edgeDir(pFace, int);
  pEdge    F_findEdge(const pFace, const pVertex, const pVertex);
  pVertex  F_edOpVt(const pFace, const pEdge);
  int      F_dirUsingEdge(pFace, pEdge);
  int      F_numRegions(pFace);
  pPList   F_regions(pFace);
  pRegion  F_region(pFace, int);
  pRegion  F_otherRegion(const pFace, const pRegion);
  int      F_inClosure(pFace, pEntity);
  void     F_chDir(pFace);
  int      F_align(pFace,pVertex,pVertex,pVertex,pVertex);

  void     F_coordP1(const pFace, double [][3]);
  void     F_coord(const pFace, double [][3]);
  bool     F_params(const pFace, double [][2]);
  double   F_area(const pFace, const double *dir=NULL);
  double   XYZ_F_area(const double [][3], const double *dir=NULL);
  double   F_areaSq(const pFace, const double *dir=NULL);
  double   XYZ_F_areaSq(const double [][3], const double *dir=NULL);
  void     F_linearParams(const pFace, const double[3], double[2]);
  void     F_center(const pFace, double [3]);
  void     F_normal(const pFace, double [3]);
  void     XYZ_F_normal(const double [3][3], double [3]);
  bool     F_volumeRatio(const pFace, double *);
  double   F_worstVolumeRatio(const std::set<pFace>);

  /********************/
  /*  Edge Operators  */
  /********************/

  void     E_info(const pEdge, std::string name="", std::ostream& out=std::cout);

  pEdge    E_exist(pVertex, pVertex);

  pGEntity E_whatIn(pEdge);
  int      E_whatInType(pEdge);
  void     E_setWhatIn(pEdge, pGEntity);

  int      E_numPoints(pEdge);
  pPoint   E_point(pEdge, int);
  pVertex  E_vertex(pEdge, int);
  pPList   E_vertices(const pEdge);
  pVertex  E_otherVertex(pEdge, pVertex);
  int      E_numFaces(pEdge);
  pPList   E_faces(const pEdge);
  pFace    E_face(pEdge, int);
  pFace    E_otherFace(pEdge, pFace);
  pFace    E_otherFace(pEdge, pFace, pRegion);
  int      E_numRegions(pEdge);
  pPList   E_regions(pEdge);
  int      E_inClosure(pEdge, pEntity);
  int      E_dir(pEdge, pVertex, pVertex);
  int      E_align(pEdge, pVertex,pVertex);

  void     E_coordP1(const pEdge, double [][3]);
  void     E_coord(const pEdge, double [][3]);
  bool     E_params(const pEdge, double [2][2]);
  double   E_length(const pEdge);
  double   E_lengthSq(const pEdge);
  double   E_linearParams(const pEdge, const pVertex);
  double   E_linearParams(const pEdge, const double[3]);
  void     E_center(const pEdge, double[3]);
  void     E_cavityCenter(const pEdge, double[3]);

  /********************/
  /* Point Operators  */
  /********************/

  double P_x(pPoint);
  double P_y(pPoint);
  double P_z(pPoint);
  void   P_setPos(pPoint, double x, double y, double z);

  void   P_setID(pPoint, int);
  int    P_id(pPoint);

  double P_param1(pPoint);
  void   P_setParam1(pPoint, double);

  /********************/
  /* Vertex Operators */
  /********************/

  void     V_info(const pVertex, std::string name="", std::ostream& out=std::cout);
 
  pGEntity V_whatIn(pVertex);
  int      V_whatInType(pVertex);
  void     V_setWhatIn(pVertex, pGEntity);

  pPoint   V_point(pVertex);
  int      V_numEdges(pVertex);
  pPList   V_edges(pVertex);
  pEdge    V_edge(pVertex, int);
  int      V_numFaces(pVertex);
  pPList   V_faces(pVertex);
  int      V_numRegions(pVertex);
  pPList   V_regions(pVertex);

  int      V_id(const pVertex);

  void     V_coord(const pVertex, double [3]);
  bool     V_params(pVertex, double *, double *);
  double   V_meanEdgeLenSq(const pVertex);
  void     V_cavityCenter(const pVertex, double[3]);

  /***************************/
  /* Entities list operators */
  /***************************/

  pPList  PList_new();
  pPList  PList_allocate( );
  void    PList_delete         (pPList);
  void    PList_deallocate     (pPList);
  void    PList_clear          (pPList);
  pPList  PList_appPListUnique (pPList, pPList source);
  pPList  PList_appUnique      (pPList, pEntity);
  pPList  PList_append         (pPList, pEntity);
  int     PList_size           (pPList); 
  pEntity PList_item           (pPList, int n);
  pEntity PList_next           (pPList, void ** restart);
  int     PList_inList         (pPList, pEntity);
  void    PList_remItem        (pPList, pEntity);
  
  /***********************/
  /* Attached data tools */
  /***********************/
  
  pMeshDataId MD_newMeshDataId(const std::string="");
  pMeshDataId MD_lookupMeshDataId(const std::string);
  void        MD_deleteMeshDataId(pMeshDataId);
  
  void EN_attachDataInt(pEntity, pMeshDataId, int);
  void EN_attachDataDbl(pEntity, pMeshDataId, double);
  void EN_attachDataPtr(pEntity, pMeshDataId, void *);
  void EN_attachDataP  (pEntity, const char *, void *);
  void EN_attachDataI  (pEntity, const char *, int);
  
  void EN_modifyDataInt(pEntity, pMeshDataId, int);
  void EN_modifyDataDbl(pEntity, pMeshDataId, double);
  void EN_modifyDataPtr(pEntity, pMeshDataId, void *);
  int  EN_modifyDataP  (pEntity, const char *, void *);
  int  EN_modifyDataI  (pEntity, const char *, int);
  
  void EN_deleteData(pEntity, pMeshDataId);
  void EN_removeData(pEntity, const char *);
  
  int  EN_getDataInt(pEntity, pMeshDataId, int *);
  int  EN_getDataDbl(pEntity, pMeshDataId, double *);
  int  EN_getDataPtr(pEntity, pMeshDataId, void **);
  void * EN_dataP(pEntity, const char *);
  int    EN_dataI(pEntity, const char *);

}

#endif 


