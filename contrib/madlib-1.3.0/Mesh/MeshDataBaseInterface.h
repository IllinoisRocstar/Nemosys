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

#include "madlib_export.h"

#include "ModelInterface.h"
// MS
#include "MeshDataBase.h"
// MS End
#include <iostream>
#include <set>

// MS
  typedef std::vector<double> v1dd;
  typedef std::vector<int> v1di;
  typedef std::vector<std::vector<double> > v2dd;
// MS End

namespace MAd {

  class MADLIB_EXPORT MDB_DataExchangerPeriodic;

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

  pMesh   MADLIB_EXPORT M_new(pGModel);
  void    MADLIB_EXPORT M_delete(pMesh);
  void    MADLIB_EXPORT M_load(pMesh, const char * filename);
  void    MADLIB_EXPORT M_clean(pMesh);

  void    MADLIB_EXPORT M_shrink(pMesh);
  void    MADLIB_EXPORT M_expand(pMesh);

  void    MADLIB_EXPORT M_info(const pMesh, std::ostream& out=std::cout);
  pGModel MADLIB_EXPORT M_model(pMesh);
  int     MADLIB_EXPORT M_dim(pMesh);

  int     MADLIB_EXPORT M_edgesMaxOrder(pMesh);
  int     MADLIB_EXPORT M_facesMaxOrder(pMesh);
  int     MADLIB_EXPORT M_regionsMaxOrder(pMesh);
  int     MADLIB_EXPORT M_maxOrder(pMesh);

  bool    MADLIB_EXPORT M_isParametric(pMesh);

  void    MADLIB_EXPORT M_writeMsh(const pMesh, const char *name, int version=2,
                                   const int * partitionTable=NULL);
  void    MADLIB_EXPORT M_writeMshPer(pMesh, const char *name,
                        MDB_DataExchangerPeriodic &, int version);
  // MS
  v2dd    MADLIB_EXPORT M_getVrtCrds(pMesh);
  v1dd    MADLIB_EXPORT M_getVrtXCrds(pMesh);
  v1dd    MADLIB_EXPORT M_getVrtYCrds(pMesh);
  v1dd    MADLIB_EXPORT M_getVrtZCrds(pMesh);
  v1di    MADLIB_EXPORT M_getVrtIds(pMesh);
  std::map<int, int> MADLIB_EXPORT M_getVrtIdxMapNew2Old(pMesh);
  std::map<int, int> MADLIB_EXPORT M_getVrtIdxMapOld2New(pMesh);
  v1di    MADLIB_EXPORT M_getConnectivities(pMesh);
  // MS End

#ifdef _HAVE_METIS_
  void    MADLIB_EXPORT M_Partition(pMesh, int nbParts, const char * filename);
#endif

  int MADLIB_EXPORT M_numRegions(pMesh);
  int MADLIB_EXPORT M_numTets(pMesh);
  int MADLIB_EXPORT M_numHexes(pMesh);
  int MADLIB_EXPORT M_numPrisms(pMesh);
  int MADLIB_EXPORT M_numFaces(pMesh);
  int MADLIB_EXPORT M_numQuads(pMesh);
  int MADLIB_EXPORT M_numTriangles(pMesh);
  int MADLIB_EXPORT M_numEdges(pMesh);
  int MADLIB_EXPORT M_numVertices(pMesh);

  RIter MADLIB_EXPORT M_regionIter(pMesh);
  FIter MADLIB_EXPORT M_faceIter(pMesh);
  EIter MADLIB_EXPORT M_edgeIter(pMesh);
  VIter MADLIB_EXPORT M_vertexIter(pMesh);

  int MADLIB_EXPORT M_numClassifiedRegions (pMesh, pGEntity);
  int MADLIB_EXPORT M_numClassifiedFaces   (pMesh, pGEntity);
  int MADLIB_EXPORT M_numClassifiedEdges   (pMesh, pGEntity);
  int MADLIB_EXPORT M_numClassifiedVertices(pMesh, pGEntity);

  RIter   MADLIB_EXPORT M_classifiedRegionIter(pMesh, pGEntity);
  FIter   MADLIB_EXPORT M_classifiedFaceIter  (pMesh, pGEntity, int closure);
  EIter   MADLIB_EXPORT M_classifiedEdgeIter  (pMesh, pGEntity, int closure);
  VIter   MADLIB_EXPORT M_classifiedVertexIter(pMesh, pGEntity, int closure);
  pVertex MADLIB_EXPORT M_classifiedVertex    (pMesh, pGVertex);

  void MADLIB_EXPORT M_classifyEntities(pMesh);

  /* mesh entity creation routines */
  pRegion MADLIB_EXPORT M_createR(pMesh, int nFace, pFace   faces[], pGEntity pg=NULL);
  pRegion MADLIB_EXPORT M_createR(pMesh, int nVert, pVertex verts[], pGEntity pg=NULL);
  pRegion MADLIB_EXPORT M_createR(pMesh, int nVert, int   vertIds[], pGEntity pg=NULL);
  pRegion MADLIB_EXPORT M_createTet(pMesh, int order, bool serendip, int vertIds[], pGEntity pg=NULL);
  pFace   MADLIB_EXPORT M_createF(pMesh, int nEdge, pEdge   edges[], pGEntity pg=NULL);
  pFace   MADLIB_EXPORT M_createF(pMesh, int nVert, pVertex verts[], pGEntity pg=NULL);
  pFace   MADLIB_EXPORT M_createF(pMesh, int nVert, int   vertIds[], pGEntity pg=NULL);
  pFace   MADLIB_EXPORT M_createTri(pMesh, int order, int vertIds[], pGEntity pg=NULL);
  pEdge   MADLIB_EXPORT M_createE(pMesh, pVertex v1, pVertex v2, pGEntity pg=NULL);
  pEdge   MADLIB_EXPORT M_createE(pMesh, int   vId1, int   vId2, pGEntity pg=NULL);
  pEdge   MADLIB_EXPORT M_createE(pMesh, pVertex v1, pVertex v2, int order=1, 
                    pVertex pts[]=NULL, pGEntity pg=NULL);
  pVertex MADLIB_EXPORT M_createV(pMesh, double x, double y, double z, int patch, pGEntity pg=NULL);
  pVertex MADLIB_EXPORT M_createV2(pMesh, double xyz[3], int patch, pGEntity pg=NULL);
  pVertex MADLIB_EXPORT M_createVP(pMesh, double x, double y, double z, 
                     double u, double v, int patch, pGEntity pg=NULL);

  /* mesh entity deletion routines */
  void   MADLIB_EXPORT M_removeRegion(pMesh, pRegion);
  void   MADLIB_EXPORT M_removeFace  (pMesh, pFace);
  void   MADLIB_EXPORT M_removeEdge  (pMesh, pEdge);
  void   MADLIB_EXPORT M_removeVertex(pMesh, pVertex);
  void   MADLIB_EXPORT P_delete(pPoint);
  pPoint MADLIB_EXPORT P_new(void);

  /* extra access to entities */
  pRegion MADLIB_EXPORT M_region(pMesh, int []);
  pFace   MADLIB_EXPORT M_face  (pMesh, int []);
  pEdge   MADLIB_EXPORT M_edge  (pMesh, int, int);
  pVertex MADLIB_EXPORT M_vertex(pMesh, int);


  /***********************/
  /* Geometric model ops */
  /***********************/
  /* deprecated (rather use ModelInterface) */

  void    MADLIB_EXPORT M_setGeomFeature    (pMesh, int, pGEntity);
  int     MADLIB_EXPORT M_numGeomFeatures   (pMesh);
  int     MADLIB_EXPORT M_geomFeatureId     (pMesh, int ith);
  pPGList MADLIB_EXPORT M_geomFeature       (pMesh, int id);
  void    MADLIB_EXPORT M_printGeomFeatures (pMesh);

  /********************/
  /* Entity Iter Ops  */
  /********************/

  pRegion MADLIB_EXPORT RIter_next  (RIter);
  void    MADLIB_EXPORT RIter_delete(RIter);
  void    MADLIB_EXPORT RIter_reset (RIter);
  pFace   MADLIB_EXPORT FIter_next  (FIter);
  void    MADLIB_EXPORT FIter_delete(FIter);
  void    MADLIB_EXPORT FIter_reset (FIter);
  pEdge   MADLIB_EXPORT EIter_next  (EIter);
  void    MADLIB_EXPORT EIter_delete(EIter);
  void    MADLIB_EXPORT EIter_reset (EIter);
  pVertex MADLIB_EXPORT VIter_next  (VIter);
  void    MADLIB_EXPORT VIter_delete(VIter);
  void    MADLIB_EXPORT VIter_reset (VIter);

  /********************/
  /* Entity Operators */
  /********************/
  
  void     MADLIB_EXPORT EN_info(const pEntity, std::string name="", std::ostream& out=std::cout);
  void     MADLIB_EXPORT EN_setID(pEntity, int id);
  int      MADLIB_EXPORT EN_id   (pEntity);
  void     MADLIB_EXPORT EN_print(pEntity);

  int      MADLIB_EXPORT EN_type(pEntity); // obsolete: returns dim, use EN_dim() instead
  int      MADLIB_EXPORT EN_dim(pEntity);
  int      MADLIB_EXPORT EN_mshTag(pEntity);
  int      MADLIB_EXPORT EN_whatInType(pEntity);
  pGEntity MADLIB_EXPORT EN_whatIn(pEntity);
  void     MADLIB_EXPORT EN_setWhatIn(pEntity,pGEntity);

  int      MADLIB_EXPORT EN_numVertices(const pEntity);
  pPList   MADLIB_EXPORT EN_vertices(const pEntity);
  int      MADLIB_EXPORT EN_inClosure(pEntity, pEntity);

  /********************/
  /* Region Operators */
  /********************/

  void     MADLIB_EXPORT R_info         (const pRegion, std::string name="", std::ostream& out=std::cout);
  void     MADLIB_EXPORT R_info_quality (const pRegion, std::ostream& out=std::cout);
  void     MADLIB_EXPORT R_info_topology(const pRegion, std::ostream& out=std::cout);

  pGRegion MADLIB_EXPORT R_whatIn(pRegion);
  int      MADLIB_EXPORT R_whatInType(pRegion);
  void     MADLIB_EXPORT R_setWhatIn(pRegion, pGEntity);

  int      MADLIB_EXPORT R_numPoints(pRegion);
  pPoint   MADLIB_EXPORT R_point(pRegion,int);
  int      MADLIB_EXPORT R_numVertices(pRegion);
  pPList   MADLIB_EXPORT R_vertices(pRegion);
  pVertex  MADLIB_EXPORT R_vertex(pRegion, int);
  pVertex  MADLIB_EXPORT R_fcOpVt(const pRegion, const pFace);
  pPList   MADLIB_EXPORT R_edges(pRegion);
  pEdge    MADLIB_EXPORT R_edge(pRegion, int);
  pEdge    MADLIB_EXPORT R_gtOppEdg(const pRegion, const pEdge);
  int      MADLIB_EXPORT R_numFaces(pRegion);
  pFace    MADLIB_EXPORT R_face(pRegion, int);
  pPList   MADLIB_EXPORT R_faces(pRegion);
  pFace    MADLIB_EXPORT R_vtOpFc(const pRegion, const pVertex);
  int      MADLIB_EXPORT R_inClosure(pRegion, pEntity);
  int      MADLIB_EXPORT R_faceDir(pRegion, int);
  int      MADLIB_EXPORT R_faceOri(pRegion, int);
  int      MADLIB_EXPORT R_dirUsingFace(pRegion,pFace);
  int      MADLIB_EXPORT R_oriUsingFace(pRegion,pFace);

  void     MADLIB_EXPORT R_coordP1(const pRegion, double [][3]);
  void     MADLIB_EXPORT R_coord(const pRegion, double [][3]);
  double   MADLIB_EXPORT R_volume(const pRegion);
  double   MADLIB_EXPORT R_XYZ_volume (const double [][3]);
  double   MADLIB_EXPORT R_circumRad(const pRegion);
  double   MADLIB_EXPORT R_inscrRad(const pRegion);
  bool     MADLIB_EXPORT R_meanRatioCube(const pRegion, double *);
  bool     MADLIB_EXPORT R_XYZ_meanRatioCube(const double [][3], double *);
  bool     MADLIB_EXPORT R_XYZ_isFlat(const double [][3]);
  void     MADLIB_EXPORT R_linearParams(const pRegion, const double[3], double[3]);
  void     MADLIB_EXPORT R_center(const pRegion, double [3]);
  void     MADLIB_EXPORT R_jacobian(const pRegion, double [3][3]);
  double   MADLIB_EXPORT R_invJacobian(const pRegion, double [3][3]);
  void     MADLIB_EXPORT R_box(const pRegion, double [3][2]);
  bool     MADLIB_EXPORT R_inBox(const pRegion, const double [3], double);
  bool     MADLIB_EXPORT R_contains(const pRegion, const double [3], double);

  /********************/
  /*  Face Operators  */
  /********************/

  void     MADLIB_EXPORT F_info(const pFace, std::string name="", std::ostream& out=std::cout);

  pFace    MADLIB_EXPORT F_exist(pEntity, pEntity, pEntity, pEntity);
  pFace    MADLIB_EXPORT F_exist(pEdge,pVertex);

  pGEntity MADLIB_EXPORT F_whatIn(pFace);
  int      MADLIB_EXPORT F_whatInType(pFace);
  void     MADLIB_EXPORT F_setWhatIn(pFace, pGEntity);

  pPoint   MADLIB_EXPORT F_point(pFace, int);
  int      MADLIB_EXPORT F_numPoints(pFace);
  int      MADLIB_EXPORT F_numPointsTot(pFace);
  int      MADLIB_EXPORT F_numVertices(pFace);
  pPList   MADLIB_EXPORT F_vertices(pFace, int dir=1);
  pVertex  MADLIB_EXPORT F_vertex(pFace, int);
  pEdge    MADLIB_EXPORT F_vtOpEd(const pFace, const pVertex);
  int      MADLIB_EXPORT F_numEdges(pFace);
  pPList   MADLIB_EXPORT F_edges(pFace);
  pEdge    MADLIB_EXPORT F_edge(pFace, int);
  int      MADLIB_EXPORT F_edgeDir(pFace, int);
  pEdge    MADLIB_EXPORT F_findEdge(const pFace, const pVertex, const pVertex);
  pVertex  MADLIB_EXPORT F_edOpVt(const pFace, const pEdge);
  int      MADLIB_EXPORT F_dirUsingEdge(pFace, pEdge);
  int      MADLIB_EXPORT F_numRegions(pFace);
  pPList   MADLIB_EXPORT F_regions(pFace);
  pRegion  MADLIB_EXPORT F_region(pFace, int);
  pRegion  MADLIB_EXPORT F_otherRegion(const pFace, const pRegion);
  int      MADLIB_EXPORT F_inClosure(pFace, pEntity);
  void     MADLIB_EXPORT F_chDir(pFace);
  int      MADLIB_EXPORT F_align(pFace,pVertex,pVertex,pVertex,pVertex);

  void     MADLIB_EXPORT F_coordP1(const pFace, double [][3]);
  void     MADLIB_EXPORT F_coord(const pFace, double [][3]);
  bool     MADLIB_EXPORT F_params(const pFace, double [][2]);
  double   MADLIB_EXPORT F_area(const pFace, const double *dir=NULL);
  double   MADLIB_EXPORT XYZ_F_area(const double [][3], const double *dir=NULL);
  double   MADLIB_EXPORT F_areaSq(const pFace, const double *dir=NULL);
  double   MADLIB_EXPORT XYZ_F_areaSq(const double [][3], const double *dir=NULL);
  void     MADLIB_EXPORT F_linearParams(const pFace, const double[3], double[2]);
  void     MADLIB_EXPORT F_center(const pFace, double [3]);
  void     MADLIB_EXPORT F_normal(const pFace, double [3]);
  void     MADLIB_EXPORT XYZ_F_normal(const double [3][3], double [3]);
  bool     MADLIB_EXPORT F_volumeRatio(const pFace, double *);
  double   MADLIB_EXPORT F_worstVolumeRatio(const std::set<pFace>);

  /********************/
  /*  Edge Operators  */
  /********************/

  void     MADLIB_EXPORT E_info(const pEdge, std::string name="", std::ostream& out=std::cout);

  pEdge    MADLIB_EXPORT E_exist(pVertex, pVertex);

  pGEntity MADLIB_EXPORT E_whatIn(pEdge);
  int      MADLIB_EXPORT E_whatInType(pEdge);
  void     MADLIB_EXPORT E_setWhatIn(pEdge, pGEntity);

  int      MADLIB_EXPORT E_numPoints(pEdge);
  pPoint   MADLIB_EXPORT E_point(pEdge, int);
  pVertex  MADLIB_EXPORT E_vertex(pEdge, int);
  pPList   MADLIB_EXPORT E_vertices(const pEdge);
  pVertex  MADLIB_EXPORT E_otherVertex(pEdge, pVertex);
  int      MADLIB_EXPORT E_numFaces(pEdge);
  pPList   MADLIB_EXPORT E_faces(const pEdge);
  pFace    MADLIB_EXPORT E_face(pEdge, int);
  pFace    MADLIB_EXPORT E_otherFace(pEdge, pFace);
  pFace    MADLIB_EXPORT E_otherFace(pEdge, pFace, pRegion);
  int      MADLIB_EXPORT E_numRegions(pEdge);
  pPList   MADLIB_EXPORT E_regions(pEdge);
  int      MADLIB_EXPORT E_inClosure(pEdge, pEntity);
  int      MADLIB_EXPORT E_dir(pEdge, pVertex, pVertex);
  int      MADLIB_EXPORT E_align(pEdge, pVertex,pVertex);

  void     MADLIB_EXPORT E_coordP1(const pEdge, double [][3]);
  void     MADLIB_EXPORT E_coord(const pEdge, double [][3]);
  bool     MADLIB_EXPORT E_params(const pEdge, double [2][2]);
  double   MADLIB_EXPORT E_length(const pEdge);
  double   MADLIB_EXPORT E_lengthSq(const pEdge);
  double   MADLIB_EXPORT E_linearParams(const pEdge, const pVertex);
  double   MADLIB_EXPORT E_linearParams(const pEdge, const double[3]);
  void     MADLIB_EXPORT E_center(const pEdge, double[3]);
  void     MADLIB_EXPORT E_cavityCenter(const pEdge, double[3]);

  /********************/
  /* Point Operators  */
  /********************/

  double MADLIB_EXPORT P_x(pPoint);
  double MADLIB_EXPORT P_y(pPoint);
  double MADLIB_EXPORT P_z(pPoint);
  void   MADLIB_EXPORT P_setPos(pPoint, double x, double y, double z);

  void   MADLIB_EXPORT P_setID(pPoint, int);
  int    MADLIB_EXPORT P_id(pPoint);

  double MADLIB_EXPORT P_param1(pPoint);
  void   MADLIB_EXPORT P_setParam1(pPoint, double);

  /********************/
  /* Vertex Operators */
  /********************/

  void     MADLIB_EXPORT V_info(const pVertex, std::string name="", std::ostream& out=std::cout);
 
  pGEntity MADLIB_EXPORT V_whatIn(pVertex);
  int      MADLIB_EXPORT V_whatInType(pVertex);
  void     MADLIB_EXPORT V_setWhatIn(pVertex, pGEntity);

  pPoint   MADLIB_EXPORT V_point(pVertex);
  int      MADLIB_EXPORT V_numEdges(pVertex);
  pPList   MADLIB_EXPORT V_edges(pVertex);
  pEdge    MADLIB_EXPORT V_edge(pVertex, int);
  int      MADLIB_EXPORT V_numFaces(pVertex);
  pPList   MADLIB_EXPORT V_faces(pVertex);
  int      MADLIB_EXPORT V_numRegions(pVertex);
  pPList   MADLIB_EXPORT V_regions(pVertex);

  int      MADLIB_EXPORT V_id(const pVertex);

  void     MADLIB_EXPORT V_coord(const pVertex, double [3]);
  bool     MADLIB_EXPORT V_params(pVertex, double *, double *);
  double   MADLIB_EXPORT V_meanEdgeLenSq(const pVertex);
  void     MADLIB_EXPORT V_cavityCenter(const pVertex, double[3]);

  /***************************/
  /* Entities list operators */
  /***************************/

  pPList  MADLIB_EXPORT PList_new();
  pPList  MADLIB_EXPORT PList_allocate( );
  void    MADLIB_EXPORT PList_delete         (pPList);
  void    MADLIB_EXPORT PList_deallocate     (pPList);
  void    MADLIB_EXPORT PList_clear          (pPList);
  pPList  MADLIB_EXPORT PList_appPListUnique (pPList, pPList source);
  pPList  MADLIB_EXPORT PList_appUnique      (pPList, pEntity);
  pPList  MADLIB_EXPORT PList_append         (pPList, pEntity);
  int     MADLIB_EXPORT PList_size           (pPList); 
  pEntity MADLIB_EXPORT PList_item           (pPList, int n);
  pEntity MADLIB_EXPORT PList_next           (pPList, void ** restart);
  int     MADLIB_EXPORT PList_inList         (pPList, pEntity);
  void    MADLIB_EXPORT PList_remItem        (pPList, pEntity);
  
  /***********************/
  /* Attached data tools */
  /***********************/
  
  pMeshDataId MADLIB_EXPORT MD_newMeshDataId(const std::string="");
  pMeshDataId MADLIB_EXPORT MD_lookupMeshDataId(const std::string);
  void        MADLIB_EXPORT MD_deleteMeshDataId(pMeshDataId);
  
  void MADLIB_EXPORT EN_attachDataInt(pEntity, pMeshDataId, int);
  void MADLIB_EXPORT EN_attachDataDbl(pEntity, pMeshDataId, double);
  void MADLIB_EXPORT EN_attachDataPtr(pEntity, pMeshDataId, void *);
  void MADLIB_EXPORT EN_attachDataP  (pEntity, const char *, void *);
  void MADLIB_EXPORT EN_attachDataI  (pEntity, const char *, int);

  void MADLIB_EXPORT EN_modifyDataInt(pEntity, pMeshDataId, int);
  void MADLIB_EXPORT EN_modifyDataDbl(pEntity, pMeshDataId, double);
  void MADLIB_EXPORT EN_modifyDataPtr(pEntity, pMeshDataId, void *);
  int  MADLIB_EXPORT EN_modifyDataP  (pEntity, const char *, void *);
  int  MADLIB_EXPORT EN_modifyDataI  (pEntity, const char *, int);

  void MADLIB_EXPORT EN_deleteData(pEntity, pMeshDataId);
  void MADLIB_EXPORT EN_removeData(pEntity, const char *);

  int  MADLIB_EXPORT EN_getDataInt(pEntity, pMeshDataId, int *);
  int  MADLIB_EXPORT EN_getDataDbl(pEntity, pMeshDataId, double *);
  int  MADLIB_EXPORT EN_getDataPtr(pEntity, pMeshDataId, void **);
  void MADLIB_EXPORT *EN_dataP(pEntity, const char *);
  int  MADLIB_EXPORT  EN_dataI(pEntity, const char *);

}

#endif 


