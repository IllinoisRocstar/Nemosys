// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Jean-Francois Remacle, Gaetan Compere, Koen Hillewaert
// -------------------------------------------------------------------

#include "MeshDataBaseInterface.h"
#include "MeshDataBaseParallelInterface.h"
#include "MeshDataBase.h"
#include "MeshDataBaseIO.h"
#include "PList.h"
#include "ParallelUtils.h"
#include "MAdDefines.h"
#include "MathUtils.h"
#include <string>
#include <map>
#ifdef PARALLEL
#include "MeshDataBaseParallelIO.h"
#endif
#include "MAdMessage.h"
#include "MAdSingleton.h"

#include <sstream>

/*! \defgroup mesh      Mesh operations */
/*! \defgroup entity    Entity operations */
/*! \defgroup region    Region operations */
/*! \defgroup face      Face operations */
/*! \defgroup edge      Edge operations */
/*! \defgroup vertex    Vertex operations */
/*! \defgroup point     Point operations */
/*! \defgroup ios       Mesh in- and output */
/*! \defgroup parallel  Communication */
/*! \defgroup internal  Internal routines */ 

typedef MAdSingleton< std::map <std::string, unsigned int> > attachableDataIds;
typedef MAdSingleton< std::map <unsigned int, std::string > > attachableDataIds_rev;

namespace MAd {

  pMeshDataId newMeshDataId(const std::string tag)
  {
    std::string tag2;
    if ( !strcmp(tag.c_str(),"") ) {
      tag2 = "X";
      bool unique = false;
      while ( !unique )
        {
          std::map<std::string, unsigned int>::iterator iter = (attachableDataIds::instance()).find(tag2);
          if(iter != (attachableDataIds::instance()).end()) {
            tag2 = tag2 + "X";
          }
          else unique = true;
        }
    }
    else tag2 = tag;

    std::map<std::string, unsigned int>::iterator iter = (attachableDataIds::instance()).find(tag2);
    if(iter != (attachableDataIds::instance()).end()) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Mesh data id with tag \'%s\' already exists",tag2.c_str());
      return (*iter).second;
    }
    if((attachableDataIds::instance()).empty())
      {
        (attachableDataIds_rev::instance())[1] = tag2;
        (attachableDataIds::instance())    [tag2] = 1;
        return 1;
      }
    else
      {
        unsigned int biggest = (*(--(attachableDataIds_rev::instance()).end())).first;
        (attachableDataIds_rev::instance())[biggest+1] = tag2;
        (attachableDataIds::instance())    [tag2] = biggest+1;
        return biggest+1;
      }
  }

  pMeshDataId MD_newMeshDataId(const std::string tag) 
  {
    return newMeshDataId(tag);
  }

  pMeshDataId MD_lookupMeshDataId(const std::string tag) {
    std::map <std::string, unsigned int>::const_iterator it = 
      (attachableDataIds::instance()).find(tag);
    if(it == (attachableDataIds::instance()).end())return newMeshDataId(tag);
    return (*it).second;
  }

  void MD_deleteMeshDataId(pMeshDataId id) {
    (attachableDataIds::instance()).erase((attachableDataIds_rev::instance())[id]);
    (attachableDataIds_rev::instance()).erase(id);
  }

  //! OBSOLETE, use EN_dim(pEntity) instead. 
  //! Returns type/dimension of element pe \ingroup entity 
  int EN_type(pEntity pe)
  {
    return pe->getDim();
  }

  //!  Returns type/dimension of element pe \ingroup entity 
  int EN_dim(pEntity pe)
  {
    return pe->getDim();
  }

  //!  Returns msh tag of element pe \ingroup entity 
  int EN_mshTag(pEntity pe)
  {
    return pe->getMshTag();
  }

  //!  Remove data with tag t from element pe \ingroup entity 

  void EN_removeData(pEntity pe, const char *t)
  {
    unsigned int itag = MD_lookupMeshDataId(t);
    pe->deleteData (itag);
  }

  //!  Add pointer data to element pe with tag \ingroup entity 
  void EN_attachDataP(pEntity pe, const char *tag, void *data)
  {
    unsigned int itag = MD_lookupMeshDataId(tag);
    EN_attachDataPtr(pe, itag, data);
  }

  //!  Remove data with tag id from element pe \ingroup entity 
  void EN_deleteData(pEntity pe, pMeshDataId id)
  {
    pe->deleteData (id);
  }

  //! Replace data with tag id with value \ingroup entity 
  void EN_modifyDataPtr(pEntity ent, pMeshDataId id, void * value)
  {
    EN_attachDataPtr(ent, id, value);
  }

  //!  Add pointer data to element pe with tag \ingroup entity 
  void EN_attachDataPtr(pEntity pe, pMeshDataId id, void * value)
  {
    mAttachableVoid *av = (mAttachableVoid *)pe->getData(id);
    if(!av)
      {
        av = new mAttachableVoid;
        pe->attachData(id,av);
      }
    av->veryNastyPointer = value;
  }

  //!  Get pointer data with tag from element \ingroup entity 
  void * EN_dataP(pEntity pe, const char *tag)
  {  
    unsigned int itag = MD_lookupMeshDataId(tag);
    mAttachableVoid *av = (mAttachableVoid *)pe->getData(itag);
    if(!av)return 0;
    return av->veryNastyPointer;
  }

  //! Replace data with tag id with value \ingroup entity 
  int EN_modifyDataP(pEntity pe, const char *tag, void * data)
  {
    EN_attachDataP(pe, tag, data);
    return 1;
  }


  //! Add integer data with tag id to element pe \ingroup entity 
  void EN_attachDataI(pEntity pe, const char *tag, int data)
  {
    pe->attachInt(MD_lookupMeshDataId(tag),data);
  }

  //! Add integer data with tag id to element pe \ingroup entity 
  void EN_attachDataInt(pEntity pe, pMeshDataId id, int data)
  {
    pe->attachInt(id,data);
  }

  //! Get integer data with tag id from element pe \ingroup entity 
  int EN_dataI(pEntity pe, const char *tag)
  {
    return pe->getAttachedInt(MD_lookupMeshDataId(tag));
  }


  //! Get floating precision data pointer with tag id from element pe \ingroup entity 
  int EN_getDataDbl(pEntity ent, pMeshDataId id, double *value)
  {   
    mAttachableDouble *ai = (mAttachableDouble *)ent->getData(id);
    if(!ai)return 0;
    *value = ai->d;
    return 1;
    //*value = ent->getAttachedDouble(id);
    //return 1;
  }

  //! Add floating precision data pointer with tag id to element pe \ingroup entity 
  void EN_attachDataDbl(pEntity ent, pMeshDataId id, double value)
  {
    ent->attachDouble(id,value);
  }


  //! Modify floating precision data pointer with tag id in element pe \ingroup entity 
  void EN_modifyDataDbl(pEntity ent, pMeshDataId id, double value)
  {
    ent->attachDouble(id,value);
  }

  //! Modify integer data with tag id in element pe \ingroup entity 
  int EN_modifyDataI(pEntity pe, const char *tag, int data)
  {
    pe->attachInt(MD_lookupMeshDataId(tag),data);
    return 1;
  }

  //! Modify integer data with tag id in element pe \ingroup entity 
  void EN_modifyDataInt(pEntity ent, pMeshDataId id, int value)
  {
    ent->attachInt(id,value);
  }


  //! Get address of data pointer with tag id in element pe \ingroup entity 
  int EN_getDataPtr(pEntity ent, pMeshDataId id, void **value)
  {
  
    mAttachableVoid *av = (mAttachableVoid *)ent->getData(id);
    if(!av)return 0;
    *value =  av->veryNastyPointer;
    return 1;
  }

  //! Get integer data with tag id from element pe \ingroup entity 
  int EN_getDataInt(pEntity ent, pMeshDataId id, int *value)
  { 
    mAttachableInt *ai = (mAttachableInt *)ent->getData(id);
    (*value) = 0;
    if(!ai)return 0;
    *value = ai->i;
    return 1;
  
    //*value = ent->getAttachedInt(id);
    //if(*value)return 1;
    //return 0;
  }

  pMesh M_new(pGModel pm)
  {
    pMesh m = new MDB_Mesh;
    m->model = pm;
    return m;
  }

  void M_delete(pMesh pm)
  {
    if (pm) { delete pm; pm=NULL; }
  }

  //! Load a mesh from a file \ingroup ios
  //!   - msh1 or msh2
  //!   - serial or parallel (format msh2 for parallel)
  //!   - periodic or non-periodic
  void M_load(pMesh pm, const char *filename)
  {
    LoadGmshMesh (pm, filename);

    pMeshDataId remoteTag =  MD_lookupMeshDataId("RemotePoint");

    V_createInfoInterface(pm,remoteTag);
    E_createInfoInterface(pm,remoteTag);
    F_createInfoInterface(pm,remoteTag);
    
  }

  //! Save a mesh \ingroup ios
  //!   - msh1 or msh2
  //!   - serial or parallel
  //!   - If a partitioning table is submitted, 
  //!     write the right partition numbers in the file
  void M_writeMsh(const pMesh mesh, const char *name, 
                  int version, const int * partitionTable)
  {
#ifdef PARALLEL
    SaveGmshMeshParallel (mesh, name, version);
#else
    SaveGmshMesh (mesh, name, version, true, partitionTable);
#endif
  }

  //! Save a periodic mesh \ingroup ios
  void M_writeMshPer(pMesh mesh, const char *name, MDB_DataExchangerPeriodic &deperiodic, int version)
  {
    SaveGmshMeshPer(mesh,name,deperiodic,version);
  }

  //! returns geometric model \ingroup mesh
  pGModel M_model(pMesh mesh)
  {
    return mesh->model;
  }

  //! reduces the mesh to its minimal datastructure  \ingroup mesh
  void M_shrink(pMesh mesh)
  {
    mesh->shrink();
  }

  //! reverts the mesh to its usable form \ingroup mesh
  void M_expand(pMesh mesh)
  {
    mesh->expand();
  }

  //! removes all deleted entities - clean up \ingroup mesh
  void M_clean(pMesh mesh)
  {
    RIter rIter = M_regionIter(mesh);
    while ( RIter_next(rIter) ) {}
    RIter_delete(rIter);

    FIter fIter = M_faceIter(mesh);
    while ( FIter_next(fIter) ) {}
    FIter_delete(fIter);
  
    EIter eIter = M_edgeIter(mesh);
    while ( EIter_next(eIter) ) {}
    EIter_delete(eIter);
  
    VIter vIter = M_vertexIter(mesh);
    while ( VIter_next(vIter) ) {}
    VIter_delete(vIter);
  }

  // -------------------------------------------------------------------
  //! Dump informations on the mesh \ingroup mesh
  void M_info(const pMesh mesh, std::ostream& out)
  {
    out << "\n";
    out << "Mesh statistics:\n";
    out << "  Number of regions : " << M_numRegions(mesh)  << "\n";
    out << "  Number of faces   : " << M_numFaces(mesh)    << "\n";
    out << "  Number of edges   : " << M_numEdges(mesh)    << "\n";
    out << "  Number of vertices: " << M_numVertices(mesh) << "\n";
    out << "\n";
  }

  // -------------------------------------------------------------------
#ifdef _HAVE_METIS_
  //! Serial function to partition a mesh and write the partitioned mesh to a file in msh2 format  \ingroup parallel
  void M_Partition(pMesh mesh, int nbParts, const char *filename)
  {
    PartitionMesh(mesh, nbParts, filename);
  }
#endif

  // -------------------------------------------------------------------
  //! Returns the dimension of the mesh \ingroup mesh
  int M_dim(pMesh pm)
  {
    if      ( M_numRegions(pm)  > 0 ) return 3;
    else if ( M_numFaces(pm)    > 0 ) return 2;
    else if ( M_numEdges(pm)    > 0 ) return 1;
    else if ( M_numVertices(pm) > 0 ) return 0;
    return -1;
  }

  // -------------------------------------------------------------------
  //! returns the maximum mapping order for mesh edges \ingroup mesh
  int M_edgeMaxOrder(pMesh mesh) 
  {
    EIter eIter = M_edgeIter(mesh);
    int o = 0;
    while (pEdge pe = EIter_next(eIter)) o = std::max(o,pe->getOrder());
    EIter_delete(eIter);
    return o;
  }
  
  //! returns the maximum mapping order for mesh face \ingroup mesh
  int M_faceMaxOrder(pMesh mesh) 
  {
    FIter fIter = M_faceIter(mesh);
    int o = 0;
    while (pFace pf = FIter_next(fIter)) o = std::max(o,pf->getOrder());
    FIter_delete(fIter);
    return o;
  }
  
  //! returns the maximum mapping order for mesh regions \ingroup mesh
  int M_regionMaxOrder(pMesh mesh) 
  {
    RIter rIter = M_regionIter(mesh);
    int o = 0;
    while (pRegion pr = RIter_next(rIter)) o = std::max(o,pr->getOrder());
    RIter_delete(rIter);
    return o;
  }

  //! returns the maximum mapping order for all elements \ingroup mesh
  int M_maxOrder(pMesh mesh) { 
    return std::max(M_edgeMaxOrder(mesh),
                    std::max(M_faceMaxOrder(mesh),
                             M_regionMaxOrder(mesh)));
  }

  // -------------------------------------------------------------------
  //! Returns true if boundary nodes have parametric coordinates \ingroup mesh
  bool M_isParametric(pMesh mesh)
  {
    return mesh->isParametric();
  }

  // -------------------------------------------------------------------
  //! returns number of regions in mesh \ingroup mesh \ingroup mesh
  int M_numRegions(pMesh pm)
  {
    return pm->nbTets + pm->nbHexes + pm->nbPrisms;
  }

  //! returns number of tetrahedra in mesh \ingroup mesh
  int M_numTets(pMesh pm)
  {
    return pm->nbTets;
  }

  //! returns number of hexahedra in mesh \ingroup mesh
  int M_numHexes(pMesh pm)
  {
    return pm->nbHexes;
  }

  //! returns number of prism in mesh \ingroup mesh
  int M_numPrisms(pMesh pm)
  {
    return pm->nbPrisms;
  }

  //! returns number of faces in mesh \ingroup mesh
  int M_numFaces(pMesh pm)
  {
    return pm->nbTriangles + pm->nbQuads;
  }

  //! returns number of triangles in mesh \ingroup mesh
  int M_numTriangles(pMesh pm)
  {
    return pm->nbTriangles;
  }

  //! returns number of quadrilaterals in mesh \ingroup mesh
  int M_numQuads(pMesh pm)
  {
    return pm->nbQuads;
  }

  //! returns number of edges in mesh \ingroup mesh
  int M_numEdges(pMesh pm)
  {
    return pm->nbEdges;
  }

  //! returns number of vertices in mesh \ingroup mesh
  int M_numVertices(pMesh pm)
  {
    return pm->nbPoints;
  }

  //! returns number of regions classified on ge \ingroup mesh
  int M_numClassifiedRegions(pMesh m,pGEntity ge) {

    return countClassifiedElements< MDB_ListT , MDB_Tet , pGEntity> (&m->tets,ge);
  }

  //! returns number of faces classified on ge \ingroup mesh
  int M_numClassifiedFaces(pMesh m,pGEntity ge) {

    return countClassifiedElements< MDB_ListF , MDB_Triangle , pGEntity> (&m->triangles,ge);
  }

  //! returns number of edges classified on ge \ingroup mesh
  int M_numClassifiedEdges(pMesh m,pGEntity ge) {

    return countClassifiedElements< MDB_ListE , MDB_Edge , pGEntity> (&m->edges,ge);
  }

  //! returns number of vertices classified on ge \ingroup mesh
  int M_numClassifiedVertices(pMesh m,pGEntity ge) {

    return countClassifiedElements< MDB_SetV , MDB_Point , pGEntity> (&m->points,ge);
  }

  //! returns region iterator over mesh \ingroup mesh
  /*! \warning not thread-safe */
  RIter M_regionIter(pMesh mesh)
  {  
    return new MDB_RegionIter (&mesh->tets,&mesh->hexes,&mesh->prisms);
  }

  //! returns face iterator over mesh \ingroup mesh
  /*! \warning not thread-safe */
  FIter M_faceIter(pMesh mesh)
  {
    return new MDB_FaceIter (&mesh->triangles,&mesh->quads);
  }

  //! returns edge iterator over mesh \ingroup mesh
  /*! \warning not thread-safe */
  EIter M_edgeIter(pMesh mesh)
  {
    return  new MDB_EIter (&mesh->edges);
  }

  //! returns vertex iterator over mesh \ingroup mesh
  /*! \warning not thread-safe */
  VIter M_vertexIter(pMesh mesh)
  {
    return new MDB_VIter (&mesh->points);
  }

  //! returns iterator for regions in \e mesh classified on model entity \e pg \ingroup mesh
  /*!
    When done with it, should be deleted  with function RIter_delete() to avoid
    memory leaks.
  */
  RIter M_classifiedRegionIter(pMesh mesh,pGEntity pg)
  {  
    return new MDB_RegionIter (&mesh->tets,&mesh->hexes,&mesh->prisms,pg);
  }

  //! returns iterator for faces in \e mesh classified on model entity \e pg \ingroup mesh
  /*!
    The argument \e c (closure) must currently be set to 0. Only the faces directly classified on \e pg will be considered. \n
    Example: if \e pg is a model region, the mesh faces classified on the model
    faces bordering \e pg will not be reachable with the present iterator. \n \n
    When done with it, should be deleted  with function FIter_delete() to avoid
    memory leaks.
  */
  FIter M_classifiedFaceIter(pMesh mesh,pGEntity pg,int c)
  {
    return new MDB_FaceIter (&mesh->triangles,&mesh->quads,pg,c);
  }

  //! returns iterator for edges in \e mesh classified on model entity \e pg \ingroup mesh
  /*!
    The argument \e c (closure) must currently be set to 0. Only the edges directly classified on \e pg will be considered. \n
    Example: if \e pg is a model region, the mesh edges classified on the model
    edges and faces bordering \e pg will not be reachable with the present iterator. \n \n
    When done with it, should be deleted  with function EIter_delete() to avoid
    memory leaks.
  */
  EIter M_classifiedEdgeIter(pMesh mesh,pGEntity pg,int c)
  {
    return  new MDB_EIter (&mesh->edges,pg,c);
  }

  //! returns iterator for vertices in \e mesh classified on model entity \e pg \ingroup mesh
  /*!
    The argument \e c (closure) must currently be set to 0. Only the vertices directly classified on \e pg will be considered. \n
    Example: if \e pg is a model region, the mesh vertices classified on the model
    vertices, edges and faces bordering \e pg will not be reachable with the present iterator. \n \n
    When done with it, should be deleted  with function VIter_delete() to avoid
    memory leaks.
  */
  VIter M_classifiedVertexIter(pMesh mesh,pGEntity pg,int c)
  {
    return new MDB_VIter (&mesh->points,pg,c);
  }

  //! returns vertex in \e mesh classified on model vertex \e pg \ingroup mesh
  /*!
    returns 0 if failed
  */
  pVertex M_classifiedVertex(pMesh mesh,pGVertex pg)
  {
    return (MDB_VIter(&mesh->points,(pGEntity) pg,0)).next();
  }

  //! Automatically classify the entities with no classification \ingroup mesh
  void M_classifyEntities(pMesh mesh)
  {
    mesh->classify_unclassified_entities();
  }

  // -------------------------------------------------------------------
  // Iterators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  pRegion RIter_next(RIter it)
  {
    return it->next();
  }
  pFace FIter_next(FIter it)
  {
    return it->next();
  }
  pEdge EIter_next(EIter it)
  {
    return it->next();
  }
  pVertex VIter_next(VIter it)
  {
    return it->next();
  }

  // -------------------------------------------------------------------
  void RIter_reset(RIter it)
  {
    it->reset();
  }
  void FIter_reset(FIter it)
  {
    it->reset();
  }
  void EIter_reset(EIter it)
  {
    it->reset();
  }
  void VIter_reset(VIter it)
  {
    it->reset();
  }

  // -------------------------------------------------------------------
  void RIter_delete(RIter it)
  {
    delete it;
  }
  void FIter_delete(FIter it)
  {
    delete it;
  }
  void EIter_delete(EIter it)
  {
    delete it;
  }
  void VIter_delete(VIter it)
  {
    delete it;
  }

  // -------------------------------------------------------------------
  //! returns the region using the vertices with these id's. \ingroup mesh
  //! Returns NULL if not found.
  pRegion M_region(pMesh m, int id[])
  {
    pFace f = M_face(m,id);
    if ( !f ) {printf("could not find face with ids: %d %d %d\n",id[0],id[1],id[2]); return 0;}

    pRegion r; pVertex v;
    for (int i=0; i<2; i++) {
      r = F_region(f,i);
      if (r) {
        v = R_fcOpVt(r,f);
        if ( v && ( V_id(v) == id[3] ) ) return r;
      }
    }

    return 0;
  }

  // -------------------------------------------------------------------
  //! returns the face using the vertices with these id's. \ingroup mesh
  //! Returns NULL if not found.
  pFace M_face(pMesh m, int id[])
  {
    pEdge e = M_edge(m,id[0],id[1]);
    if ( !e ) {printf("could not find edge with ids: %d %d\n",id[0],id[1]); return 0;}

    pPList eF = E_faces(e);
    void * tmp = NULL;
    pFace f;
    pVertex v;
    while ( ( f = (pFace)PList_next(eF,&tmp) ) ) {
      v = F_edOpVt(f,e);
      if ( v && ( V_id(v) == id[2] ) ) {
        PList_delete(eF);
        return f;
      }
    }

    PList_delete(eF);
    return 0;
  }

  // -------------------------------------------------------------------
  //! returns the edge between the vertices with these id's. \ingroup mesh
  //! Returns NULL if not found.
  pEdge M_edge(pMesh m, int id0, int id1)
  {
    return m->find_edge(id0,id1);
  }

  // -------------------------------------------------------------------
  //! returns the vertex with this id. Returns NULL if not found. \ingroup mesh
  pVertex M_vertex(pMesh m, int id)
  {
    return m->find_point(id);
  }

  // -------------------------------------------------------------------
  // Region operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  //! Dumps all informations about the region. \ingroup region
  void R_info(const pRegion region, std::string name, std::ostream& out)
  {
    out << "\nRegion \'" << name << "\' (" << region 
        << ", Id: "<<(int)EN_id((pEntity)region)
        <<") informations:\n";

    out << "  Classification:  ";
    pGEntity pGE = (pGEntity)R_whatIn(region);
    if (!pGE) {
      out << "NULL";
    }
    else {
      out << "GEntity: " << pGE 
          << ", dim: " << GEN_type(pGE)
          << ", tag: " << GEN_tag(pGE);
    }
    out << "\n";

    R_info_quality(region, out);
    R_info_topology(region, out);

    out << "\n";
  }

  // -------------------------------------------------------------------
  //! Dumps quality informations about the region. \ingroup region
  void R_info_quality(const pRegion region, std::ostream& out)
  {
    out << "  Quality informations:\n";
    double volume = R_volume(region);
    out << "    Volume      : " << volume  << "\n";
    if (volume < 0.) out << "  *** Negative volume ***\n";
    else {
      out << "    r/R ratio : " << R_inscrRad(region)/R_circumRad(region)<< "\n";
      double meanRatio3;
      R_meanRatioCube(region,&meanRatio3);
      out << "    Cubic mean ratio : " << meanRatio3 << "\n";
    }
  }

  // -------------------------------------------------------------------
  //! Dumps topology informations about the region. \ingroup region
  void R_info_topology(const pRegion region, std::ostream& out)
  {
    out << "  Topology informations:\n";
    out << "\n    Faces       (ids):   type  geom ent.      tag   orient\n";
    for (int iF=0; iF<R_numFaces(region); iF++) {
      pFace face = R_face(region,iF);
      pGEntity pGE = F_whatIn(face);
      int gTag = pGE ? GEN_tag(pGE) : -1;
      int gDim = pGE ? GEN_type(pGE) : -1;
      out << " " << face << " (" << EN_id((pEntity)face) << ") " 
          << gDim << "   " << pGE << " " << gTag << " " 
          << R_faceDir(region,iF) << "\n";
    }
    out << "\n    Edges       (ids):   type  model ent.      tag   orient\n";
    pPList rEdges = R_edges(region);
    void * temp = NULL;
    while( pEdge edge = (pEdge)PList_next(rEdges,&temp) ) {
      pGEntity pGE = E_whatIn(edge);
      int gTag = pGE ? GEN_tag(pGE) : -1;
      int gDim = pGE ? GEN_type(pGE) : -1;
      out << " " << edge << " (" << EN_id((pEntity)edge) << ") " 
          << gDim << "   " << pGE << " " << gTag << "\n";
    }
    PList_delete(rEdges);

    out << "\nVertices    (ids):   type  model ent.      tag               coordinates\n";
    pPList rVerts = R_vertices(region);
    temp = NULL;
    while( pVertex pV = (pVertex)PList_next(rVerts,&temp) ) {
      double xyz[3];
      V_coord(pV,xyz);
      pGEntity pGE = V_whatIn(pV);
      int gTag = pGE ? GEN_tag(pGE) : -1;
      int gDim = pGE ? GEN_type(pGE) : -1;
      out << " " << pV << " (" << EN_id((pEntity)pV) << ") " 
          << gDim << "   " << pGE << " " << gTag << " " 
          << xyz[0] << " " << xyz[1] << " " << xyz[2] << "\n";
    }
    PList_delete(rVerts);
  }

  // -------------------------------------------------------------------
  //! Returns number of faces for region "r" \ingroup region
  int R_numFaces(pRegion pr)
  {
    return pr->getNbFace ();
  }

  // -------------------------------------------------------------------
  //! Returns n-th face for region "r" \ingroup region
  pFace R_face(pRegion pr, int n)
  {
    return pr->getFace(n);
  }

  // -------------------------------------------------------------------
  //! Returns number of edges for region "r" \ingroup region
  int R_numEdges(pRegion pr)
  {
    return pr->getNbEdge();
  }

  // -------------------------------------------------------------------
  //! Classify region "r" on geometric entity "ge" \ingroup region
  void R_setWhatIn(pRegion region, pGEntity what)
  {
    region->g = what;
    for (int iF=0; iF<region->getNbFace(); iF++) {
      pFace face = region->getFace(iF);
      if ( GEN_type(face->g) > GEN_type(what) ) F_setWhatIn(face,what);
    }
  }

  // -------------------------------------------------------------------
  //! Returns n-th edge of region "r"
  pEdge R_edge(pRegion pr, int n)
  { 
    return pr->getEdge(n);
  }

  // -------------------------------------------------------------------
  //! Returns number of principal vertices in region "r" \ingroup region
  int R_numVertices(pRegion pr)
  {
    return pr->getNbVertex();
  }

  // -------------------------------------------------------------------
  //! Returns nth vertex in region "r" \ingroup region
  pVertex R_vertex(pRegion pr, int n)
  {
    return pr->getVertex (n);
  }

  // -------------------------------------------------------------------
  //! Returns 1 if the face normal points outwards of region pr, 0 otherwise \ingroup region
  int R_faceDir (pRegion pr, int n)
  {
    return pr->getFaceDir(n);
  }

  // -------------------------------------------------------------------
  //! Determine face orientation with respect to the template of the region "r" \ingroup region
  /*!
    \warning Only implemented for tetrahedra
    Returns s*(n+1), if the face is part of the closure of the region, 0 otherwise \n
    Here
    \li n is the number of times we need to rotate the face to correspond 
    with the principal vertex of the face in the element
    \li s is the sign of the normal with respect to that of the template
  */
  int R_faceOri(pRegion pr,int n) 
  {
    return pr->getFaceOrientation(n);
  }

  // -------------------------------------------------------------------
  //! Return a list of ordered edges \ingroup region
  pPList R_edges(pRegion pr)
  {
    pPList pl = PList_new();
    for (int i=0;i<pr->getNbEdge();i++)
      PList_append(pl,R_edge(pr,i));
    return pl;
  }

  // -------------------------------------------------------------------
  //! Returns a list of ordered faces \ingroup region
  pPList R_faces(pRegion pr)
  {
    pPList pl = PList_new();
    for (int i=0;i<pr->getNbFace();i++)
      PList_append(pl,R_face(pr,i));
    return pl;
  }

  // -------------------------------------------------------------------
  //! Returns a list of ordered principal vertices \ingroup region
  pPList R_vertices(pRegion pr)
  {
    pPList pl = PList_new();
    for (int i=0;i<pr->getNbVertex();i++)
      PList_append(pl,R_vertex(pr,i));
    return pl;
  }

  // -------------------------------------------------------------------
  //! Verify whether or not entity ent is a principal vertex, edge or face of region "r" \ingroup region
  int R_inClosure(pRegion pr, pEntity ent)
  {
    int dim = EN_type (ent);
    if (dim == 2)
      {
        for (int i=0;i<pr->getNbFace();i++)
          {
            if (R_face(pr,i)== ent) return 1;
          }
        return 0;
      }
    else if (dim ==1)
      {
        for (int i=0;i<pr->getNbEdge();i++)
          {
            if (R_edge(pr,i)== ent) return 1;
          }
        return 0;
      }
    else if (dim ==0)
      {
        for (int i=0;i<pr->getNbVertex();i++)
          {
            if (R_vertex(pr,i) == ent) return 1;
          }
        return 0;
      }
    else throw;
  }

  // -------------------------------------------------------------------
  // !return 1 if the face direction points to outside of the tet
  // !return 0                                 inside
  int R_dirUsingFace(pRegion pr, pFace face)
  {
    for (int i=0;i<pr->getNbFace();i++)
      if (pr->getFace(i) == face)
        return R_faceDir (pr,i);
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns the orientation of the face pf in region pr \ingroup region
  /*! Returns -1 if the face is not in the closure of the region */
  int R_oriUsingFace(pRegion pr, pFace face)
  {
    for (int i=0;i<pr->getNbFace();i++)
      if (pr->getFace(i) == face)
        return R_faceOri (pr,i);
    return -1;
  }

  // -------------------------------------------------------------------
  //! Returns n-th higher order point inside of the region, excluding those in its closure \ingroup region
  pPoint R_point(pRegion pr, int n)
  {
    return pr->getHighOrderPoint (n);
  }

  // -------------------------------------------------------------------
  //! Returns the number of higher order points of the region, excluding those in its closure \ingroup region
  int R_numPoints(pRegion pr)
  {
    return pr->getNbHighOrderPoints();
  }

  // -------------------------------------------------------------------
  //! Returns 1 if the direction of the edge in pf follows the template, 0 otherwise \ingroup face
  int F_dirUsingEdge(pFace pf, pEdge edge)
  {
    const int nbEdge = pf->getNbEdges ();
    for (int i=0;i<nbEdge;i++)
      if (pf->getEdge(i) == edge)
        return F_edgeDir (pf,i);
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns the geometrical entity on which the region "r" is classified \ingroup region
  pGRegion R_whatIn(pRegion pe)
  {
    return (pGRegion) pe->g;
  }

  // -------------------------------------------------------------------
  //! Returns the type of the geometrical entity on which region "r" is classified \ingroup region
  int R_whatInType(pRegion e) { return EN_whatInType(e);}

  // -------------------------------------------------------------------
  //! Get coordinates of vertices of the region not including high-order points. \ingroup region
  void R_coordP1(const pRegion region, double xyz[][3])
  {
    int iNode = 0;

    // Summits of the region (first order nodes)
    pPList rVerts = R_vertices(region);
    void * temp = NULL;
    while ( pVertex pV = (pVertex)PList_next(rVerts,&temp) ) {
      V_coord(pV,xyz[iNode]);
      iNode++;
    }
    PList_delete(rVerts);
  }

  // -------------------------------------------------------------------
  //! Get coordinates of vertices and high-order points of the region. \ingroup region
  void R_coord(const pRegion region, double xyz[][3])
  {
    // Summits of the region (first order nodes)
    R_coordP1(region,xyz);
    int iNode = region->getNbVertex();

    // points on edges (higher order nodes)
    pPList rEdges = R_edges(region);
    void * temp = NULL;
    while ( pEdge edge = (pEdge) PList_next(rEdges, &temp)) {
      int nEPts = E_numPoints(edge);
      for (int iEP=0; iEP<nEPts; iEP++) {
        pPoint pP = E_point(edge,iEP);
        xyz[iNode][0] = P_x(pP); xyz[iNode][1] = P_y(pP); xyz[iNode][2] = P_z(pP);
        iNode++;
      }
    }
    PList_delete(rEdges);

    // points on face (higher order nodes)
    pPList rFaces = R_faces(region);
    temp = NULL;
    while ( pFace face = (pFace) PList_next(rFaces, &temp)) {
      int nFPts = F_numPoints(face);
      for (int iFP=0; iFP<nFPts; iFP++) {
        pPoint pP = F_point(face,iFP);
        xyz[iNode][0] = P_x(pP); xyz[iNode][1] = P_y(pP); xyz[iNode][2] = P_z(pP);
        iNode++;
      }
    }
    PList_delete(rFaces);

    // points on region (higher order nodes)
    int nRPts = R_numPoints(region);
    for (int iRP=0; iRP<nRPts; iRP++) {
      pPoint pP = R_point(region,iRP);
      xyz[iNode][0] = P_x(pP); xyz[iNode][1] = P_y(pP); xyz[iNode][2] = P_z(pP);
      iNode++;
    }
  }

  // -------------------------------------------------------------------
  //! Returns the physical volume of the region. \ingroup region
  double R_volume(pRegion region)
  {
    double xyz[12][3];
    R_coordP1(region,xyz);
    return R_XYZ_volume(xyz);
  }

  // -------------------------------------------------------------------
  //! Returns the physical volume of the region with coordinates xyz. \ingroup region 
  double R_XYZ_volume (const double xyz[][3])
  {
    double e01[3], e02[3], e03[3];
    diffVec(xyz[1],xyz[0],e01);
    diffVec(xyz[2],xyz[0],e02);
    diffVec(xyz[3],xyz[0],e03);
    double nor012[3];
    crossProd(e01,e02,nor012);
    return ( dotProd(nor012,e03) * MAdSIXTH );
  }

  // -------------------------------------------------------------------
  //! Returns the circumradius of the region \ingroup region
  double R_circumRad(const pRegion region)
  {
    double xyz[4][3];
    R_coordP1(region,xyz);

    double edges0x[3][3];
    for (int i=0; i<3; i++) diffVec(xyz[i+1],xyz[0],edges0x[i]);
    double detEdgesInv = 1. / ( detMat(edges0x) );

    double edgeMids[3][3];
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        edgeMids[i][j] = (xyz[0][j] + xyz[i+1][j]) * 0.5;
      }
    }

    double tmpVec[3];
    for (int i=0; i<3; i++) tmpVec[i] = dotProd(edges0x[i],edgeMids[i]);
  
    double center[3];
    double tmpMat[3][3];
    for (int i=0; i<3; i++) {
      tmpMat[i][0] = tmpVec[i];
      tmpMat[i][1] = edges0x[i][1];
      tmpMat[i][2] = edges0x[i][2];
    }
    center[0] = detMat(tmpMat);

    for (int i=0; i<3; i++) {
      tmpMat[i][0] = edges0x[i][0];
      tmpMat[i][1] = tmpVec[i];
      tmpMat[i][2] = edges0x[i][2];
    }
    center[1] = detMat(tmpMat);

    for (int i=0; i<3; i++) {
      tmpMat[i][0] = edges0x[i][0];
      tmpMat[i][1] = edges0x[i][1];
      tmpMat[i][2] = tmpVec[i];
    }
    center[2] = detMat(tmpMat);
    
    for (int i=0; i<3; i++) { center[i] *= detEdgesInv; }
  
    double tmpVec2[3];
    diffVec(center,xyz[0],tmpVec2);

    return sqrt( tmpVec2[0]*tmpVec2[0] + tmpVec2[1]*tmpVec2[1] + tmpVec2[2]*tmpVec2[2]);
  }

  // -------------------------------------------------------------------
  //! Returns the inscribed radius of the region \ingroup region
  double R_inscrRad(const pRegion region)
  {
    double xyz[4][3];
    R_coordP1(region,xyz);

    double edges0x[5][3];
    for (int i=0; i<3; i++) {
      diffVec(xyz[i+1],xyz[0],edges0x[i]);
    } 
    diffVec(xyz[1],xyz[2],edges0x[3]);
    diffVec(xyz[2],xyz[3],edges0x[4]);

    double tmpVec[3];
    double A = 0.;

    crossProd(edges0x[0],edges0x[1],tmpVec);
    A += sqrt( dotProd(tmpVec,tmpVec) ) * 0.5;
    double V = dotProd(tmpVec,edges0x[2]) * MAdSIXTH;
    crossProd(edges0x[1],edges0x[2],tmpVec);
    A += sqrt( dotProd(tmpVec,tmpVec) ) * 0.5;
    crossProd(edges0x[2],edges0x[0],tmpVec);
    A += sqrt( dotProd(tmpVec,tmpVec) ) * 0.5;
    crossProd(edges0x[3],edges0x[4],tmpVec);
    A += sqrt( dotProd(tmpVec,tmpVec) ) * 0.5;

    return ( 3. * V ) / A;
  }

  // -------------------------------------------------------------------
  //! Computes the cubic mean ratio of a region. Returns 0 if negative volume. \ingroup region
  bool R_meanRatioCube(const pRegion region, double * mrc)
  {
    double xyz[4][3];
    R_coordP1(region,xyz);
    return R_XYZ_meanRatioCube(xyz,mrc) ;
  }

  // -------------------------------------------------------------------
  //! Computes the cubic mean ratio of a region with coordinates xyz. Returns false if negative volume. \ingroup region
  bool R_XYZ_meanRatioCube(const double xyz[][3], double * mrc)
  {
    double edges[6][3];
    diffVec(xyz[0],xyz[1],edges[0]);
    diffVec(xyz[0],xyz[2],edges[1]);
    diffVec(xyz[0],xyz[3],edges[2]);
    diffVec(xyz[1],xyz[2],edges[3]);
    diffVec(xyz[1],xyz[3],edges[4]);
    diffVec(xyz[2],xyz[3],edges[5]);

    // compute the sum of edges length square
    double lSq = 0.;
    for (int iE=0; iE<6; iE++) {
      lSq += dotProd(edges[iE],edges[iE]);
    }

    // compute volume
    double vol = R_XYZ_volume(xyz);
    if ( vol < 0. ) {
      *mrc = 0.;
      return false;
    }

    // compute cubic mean ratio
    *mrc = 15552. * vol * vol / ( lSq * lSq * lSq );
    return true;
  }

  // -------------------------------------------------------------------
  //! Check if the region with coordinates xyz \ingroup region
  //! is nearly flat or with a negative volume.
  bool R_XYZ_isFlat(const double xyz[][3])
  {
    for(int i=0; i<4; i++) {
      int ind[3];
      switch (i) {
      case 0:
        ind[0] = 1;
        ind[1] = 3;
        ind[2] = 2;
        break;
      case 1:
        ind[0] = 0;
        ind[1] = 2;
        ind[2] = 3;
        break;
      case 2:
        ind[0] = 3;
        ind[1] = 1;
        ind[2] = 0;
        break;
      case 3:
        ind[0] = 0;
        ind[1] = 1;
        ind[2] = 2;
        break;
      }

      double v01[3], v02[3];
      diffVec(xyz[ind[1]],xyz[ind[0]],v01);
      diffVec(xyz[ind[2]],xyz[ind[0]],v02);

      double normal[3];
      crossProd(v01,v02,normal);
      double ASq = dotProd(normal,normal);

      double v0X[3];
      diffVec(xyz[i],xyz[ind[0]],v0X);
      double distA = dotProd(v0X,normal);
    
      if( distA <= 0. || distA*distA < MAdTOL*MAdTOL*ASq) return true;
    }

    return false;
  }

  // -------------------------------------------------------------------
  //! Returns the vertex of the region opposite to the face. \ingroup region
  pVertex R_fcOpVt(const pRegion region, const pFace face)
  {
    pVertex opp = NULL;

    pPList fVerts = F_vertices(face,1);
    pPList rVerts = R_vertices(region);

    void * tempR = NULL;
    while( pVertex pRV = (pVertex)PList_next(rVerts,&tempR) ) {
      bool foundVertInFace = false;
      void * tempF = NULL;
      while( pVertex pFV = (pVertex)PList_next(fVerts,&tempF) ) {
        if ( pFV == pRV ) {
          foundVertInFace = true;
          break;
        }
      }
      if ( !foundVertInFace ) {
        opp = pRV;
        break;
      }
    }
    PList_delete(rVerts);
    PList_delete(fVerts);

    return opp;
  }

  // -------------------------------------------------------------------
  //! Returns the edge of the region opposite to given edge. \ingroup region
  pEdge R_gtOppEdg(const pRegion region, const pEdge edge)
  {
    assert( region->getNbEdge() == 6 );

    pVertex v0 = E_vertex(edge,0);
    pVertex v1 = E_vertex(edge,1);

    for (int iE=0; iE<6; iE++) {

      pEdge rEdge = R_edge(region,iE);

      if ( rEdge == edge ) continue;

      pVertex rEV0 = E_vertex(rEdge,0);
      if ( rEV0 == v0 || rEV0 == v1 ) continue;
      pVertex rEV1 = E_vertex(rEdge,1);
      if ( rEV1 == v0 || rEV1 == v1 ) continue;

      return rEdge;
    }

    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "could not find opposite edge");
    throw;
    return NULL;
  }

  // -------------------------------------------------------------------
  //! Returns the face of the region opposite to the vertex. \ingroup region
  pFace R_vtOpFc(const pRegion region, const pVertex vertex)
  {
    assert( region->getNbFace() == 4 );
  
    for (int iF=0; iF<R_numFaces(region); iF++ ) {
      pFace face = R_face(region,iF);
      pPList fVerts = F_vertices(face,1);
      if ( PList_inList(fVerts,vertex) ) {
        PList_delete(fVerts);
        continue;
      }
      PList_delete(fVerts);
      return face;
    }

    return NULL;
  }
  
  // -------------------------------------------------------------------
  //! Finds the coordinates in the parent element 
  void R_linearParams(const pRegion pR, const double xyz[3], 
                      double res[3])
  {
    double rxyz[4][3];
    R_coordP1(pR,rxyz);
  
    double mat[3][3];
    mat[0][0] = rxyz[1][0] - rxyz[0][0];
    mat[0][1] = rxyz[1][1] - rxyz[0][1];
    mat[0][2] = rxyz[1][2] - rxyz[0][2];
    mat[1][0] = rxyz[2][0] - rxyz[0][0];
    mat[1][1] = rxyz[2][1] - rxyz[0][1];
    mat[1][2] = rxyz[2][2] - rxyz[0][2];
    mat[2][0] = rxyz[3][0] - rxyz[0][0];
    mat[2][1] = rxyz[3][1] - rxyz[0][1];
    mat[2][2] = rxyz[3][2] - rxyz[0][2];
    
    double invMat[3][3];
    inverseMat(mat,invMat);
    
    double vec[3];
    diffVec(xyz,rxyz[0],vec);

    vecMat(vec,invMat,res);
  }

  // -------------------------------------------------------------------
  //! Returns the center of a linear region. \ingroup region
  void R_center(const pRegion region, double center[3])
  {
    double rxyz[4][3];
    R_coordP1(region,rxyz);
    center[0] = 0.25 * ( rxyz[0][0] + rxyz[1][0] + rxyz[2][0] + rxyz[3][0] );
    center[1] = 0.25 * ( rxyz[0][1] + rxyz[1][1] + rxyz[2][1] + rxyz[3][1] );
    center[2] = 0.25 * ( rxyz[0][2] + rxyz[1][2] + rxyz[2][2] + rxyz[3][2] );
  }

  // -------------------------------------------------------------------
  // Get the Jacobian of a tetrahedron
  void R_jacobian(const pRegion region, double jac[3][3])
  {
    double xyz[4][3];
    R_coordP1(region,xyz);
    for (int i=0; i<3; i++) {
      jac[i][0] = xyz[1][i] - xyz[0][i];
      jac[i][1] = xyz[2][i] - xyz[0][i];
      jac[i][2] = xyz[3][i] - xyz[0][i];
    }
  }

  // -------------------------------------------------------------------
  // Get the inverse of the Jacobian of a tetrahedron
  // returns the determinant of the Jacobian
  double R_invJacobian(const pRegion region, double ijac[3][3])
  {
    double jac[3][3];
    R_jacobian(region,jac);
    return inverseMat(jac,ijac);
  }

  // -------------------------------------------------------------------
  //! returns the bounding box (xyz-oriented) of the region \ingroup region
  //! warning: only implemented for tets
  void R_box(const pRegion pr, double box[3][2])
  {
    double xyz[4][3];
    R_coordP1(pr,xyz);
    
    for (int i=0; i<3; i++) {
      box[i][0] = MAdBIG;
      box[i][1] = -MAdBIG;
      for (int j=0; j<4; j++) {
        box[i][0] = std::min(box[i][0],xyz[j][i]);
        box[i][1] = std::max(box[i][1],xyz[j][i]);
      }
    }
  }

  // -------------------------------------------------------------------
  //! returns true if the point is in the bounding box (xyz-oriented)
  //! of the region plus a tolerance \ingroup region
  //! warning: only implemented for tets
  bool R_inBox(const pRegion pr, const double xyz[3], double tol)
  {
    double box[3][2];
    R_box(pr,box);

    for (int i=0; i<3; i++) {
      if ( xyz[i] < box[i][0] - tol || xyz[i] > box[i][1] + tol ) return false;
    }
    
    return true;
  }

  // -------------------------------------------------------------------
  //! returns true if the point is in the region plus a tolerance \ingroup region
  //! warning: only implemented for tets
  bool R_contains(const pRegion pr, const double xyz[3], double tol)
  {
    double tmp[4];
    R_linearParams(pr,xyz,tmp);
    tmp[3] = 1. - tmp[0] - tmp[1] - tmp[2];

    for (int i=0; i<4; i++) {
      if ( tmp[i] < -tol || tmp[i] > 1.+tol ) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Face operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  //! Dumps informations about the face. \ingroup face
  void F_info(const pFace face, std::string name, std::ostream& out)
  {
    out << "\n";
    out << "Face \'" << name << "\' " << face
        << ", id: " << EN_id((pEntity)face) << "\n";
    out << "  Classified on " << face->g 
        << " with dim " << GEN_type(face->g) 
        << " and tag "  << GEN_tag(face->g) << "\n";

    out << "\n--- Regions:\n";
    for (int iR=0; iR<2; iR++) {
      pRegion region = face->getRegion(iR);
      if (!region) {
        out << "Region "<<iR<<": "<<region<<"\n";
      }
      else {
        out << "Region "<<iR<<" (" << region 
            << ", Id: "<<(int)EN_id((pEntity)region)
            <<") informations:\n";
        out << "  Classification:  ";
        pGEntity pGE = (pGEntity)R_whatIn(region);
        if (!pGE) {
          out << "NULL";
        }
        else {
          out << "GEntity: " << pGE 
              << ", dim: " << GEN_type(pGE)
              << ", tag: " << GEN_tag(pGE);
        }
        out << "  Vertex opposite to face: \n";
        V_info(R_fcOpVt(region,face),"",out);
        out << "\n";
      }
    }

    out << "\n--- Edges:\n";
    for (int iE=0; iE<face->getNbEdges(); iE++) {
      pEdge edge = face->getEdge(iE);
      out << "Edge " << iE << ":\n";
      E_info(edge, "", out);
    }

//       out << "Edge "<<iE<<" ("<< edge <<"), id " 
//           << EN_id((pEntity)edge) <<":\n";
//       pGEntity pGE = EN_whatIn((pEntity)edge);
//       out << "  - Classification: dim: " << GEN_type(pGE) << ", tag: " << GEN_tag(pGE) << "\n";
//       pVertex pv0 = E_vertex(edge,0);
//       pVertex pv1 = E_vertex(edge,1);
//       pGEntity pGE0 = EN_whatIn((pEntity)pv0);
//       pGEntity pGE1 = EN_whatIn((pEntity)pv1);
//       out << "  - Vertices classifications:\n"
//           << "      - V0: dim: " << GEN_type(pGE0) << ", tag: " << GEN_tag(pGE0) << "\n"
//           << "      - V1: dim: " << GEN_type(pGE1) << ", tag: " << GEN_tag(pGE1) << "\n";
//       out << "\n";
//     }
    out << "\n";
  }

  // -------------------------------------------------------------------
  //! Returns number of edges in pf \ingroup face
  int F_numEdges(pFace pf)
  {
    return pf->getNbEdges();
  }

  // -------------------------------------------------------------------
  //! Returns n-th edge in pf \ingroup face
  pEdge F_edge(pFace pf, int n)
  {
    return pf->getEdge(n);
  }

  // -------------------------------------------------------------------
  //! Returns the edge used by the face and using the two vertices
  //! Returns NULL if not found \ingroup face
  pEdge F_findEdge(const pFace pf, const pVertex v0, const pVertex v1)
  {
    return pf->find_edge(v0,v1);
  }

  // -------------------------------------------------------------------
  //! Returns number of vertices in pf \ingroup face
  int F_numVertices(pFace pf)
  {
    return pf->getNbNodes();
  }

  // -------------------------------------------------------------------
  //! Returns n-th vertex in pf \ingroup face
  pVertex F_vertex(pFace pf, int n)
  {
    return pf->getNode(n);
  }

  // -------------------------------------------------------------------
  //! Returns 1 if the direction of the n-th edge in pf follows the template, 0 otherwise \ingroup face
  int F_edgeDir (pFace pf, int n)
  {
    pEdge e = F_edge (pf,n);
    pVertex ve1 = E_vertex (e,0);
    pVertex ve2 = E_vertex (e,1);
    pVertex vf1 = F_vertex (pf,n);
    pVertex vf2 = F_vertex (pf,(n+1)%F_numVertices(pf));
    if (ve1 == vf1 && ve2 == vf2) return 1;
    if (ve2 == vf1 && ve1 == vf2) return 0;
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns number of regions referring to pf \ingroup face
  int F_numRegions(pFace pf)
  {
    return pf->getNbRegions();
  }

  // -------------------------------------------------------------------
  //! Returns the ordered list of vertices of pf \ingroup face
  //! if dir \< 0 order is inverted, else the template is followed
  pPList F_vertices(pFace pf, int dir)
  {
    pPList pl = PList_new();
    int n = F_numVertices(pf);
    if (dir<=0)
      {
        for (int i=n-1;i>=0;i--)
          PList_append(pl,F_vertex(pf,i));
      }
    else
      {
        for (int i=0;i<n;i++)
          PList_append(pl,F_vertex(pf,i));
      }
    return pl;
  }

  // -------------------------------------------------------------------
  //! Returns the ordered list of edges composing the closure of pf \ingroup face
  pPList F_edges(pFace pf)
  {
    pPList pl = PList_new();
    int n = F_numVertices(pf);
    for (int i=0;i<n;i++)
      PList_append(pl,F_edge(pf,i));
    return pl;
  }

  // -------------------------------------------------------------------
  //! Returns the n-th region attached to pf \ingroup face
  pRegion F_region(pFace pf, int n)
  {
    return pf->getRegion(n);
  }

  // -------------------------------------------------------------------
  //! Returns the list of regions attached to pf \ingroup face
  pPList F_regions(pFace pf)
  {
    pPList pl = PList_new();
    for (int i=0;i<pf->getNbRegions();++i)
      PList_append(pl,F_region(pf,i));  
    return pl;
  }

  // -------------------------------------------------------------------
  //! Returns the region around 'pf' which is not 'pr' or NULL if none \ingroup face
  pRegion F_otherRegion(const pFace pf, const pRegion pr)
  {
    if ( pf->getRegion(0) == pr ) return pf->getRegion(1);
    if ( pf->getRegion(1) == pr ) return pf->getRegion(0);
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Face not in the closure of region");
    return NULL;
  }

  // -------------------------------------------------------------------
  //! Returns 1 if ent is in the closure of pf, else 0 \ingroup face
  int F_inClosure(pFace f, pEntity ent)
  {
    int dim = EN_type (ent);

    if (dim ==1)
      {
        for (int i=0;i<f->getNbEdges();i++)
          {
            if (F_edge(f,i)== ent) return 1;
          }
        return 0;
      }
    else if (dim ==0)
      {
        for (int i=0;i<f->getNbEdges();i++)
          {
            if (F_vertex(f,i) == ent) return 1;
          }
        return 0;
      }
    else throw;
  }

  // -------------------------------------------------------------------
  //! Returns the geometric entity to which pf is attached \ingroup face
  pGEntity F_whatIn(pFace pf)
  {
    return pf->g;
  }

  // -------------------------------------------------------------------
  //! Returns the type/dimension of the geometric entity to which pf is attached \ingroup face
  int F_whatInType(pFace e)   { return EN_whatInType(e);}

  // -------------------------------------------------------------------
  //! Classify pf on geometric entity 'what' \ingroup face
  void F_setWhatIn(pFace face, pGEntity what)
  {
    face->g = what;
    for (int iE=0; iE<face->getNbEdges(); iE++) {
      pEdge edge = face->getEdge(iE);
      if ( GEN_type(edge->g) > GEN_type(what) ) E_setWhatIn(edge,what);
    }
  }

  // -------------------------------------------------------------------
  //! Inverts the direction of face pf by switching edge 1 and 2 \ingroup face
  /*! \warning only implemented for triangles */
  void F_chDir(pFace pf) 
  {
    assert (pf->getNbEdges() == 3);
    {
      MDB_Triangle *t = (MDB_Triangle*)pf;
      pEdge temp = t->e2;
      t->e2 = t->e1;
      t->e1 = temp;
    }
  }

  // -------------------------------------------------------------------
  //! Returns the n-th higher-order point attached to pf, excluding the closure \ingroup face
  pPoint F_point(pFace f, int n)
  {
    return f->getHighOrderPoint (n);
  }

  // -------------------------------------------------------------------
  //! Returns the number of higher-order points attached to pf, excluding the closure \ingroup face
  int F_numPoints(pFace f)
  {
    return f->getNbHighOrderPoints();
  }

  // -------------------------------------------------------------------
  //! Returns the number of points attached to pf, including summits and points in closure \ingroup face
  int F_numPointsTot(pFace f)
  {
    int numPt = f->getNbNodes() + f->getNbHighOrderPoints();

    pPList list = F_edges(f);
    void * temp = NULL;
    while ( pEdge edge = (pEdge) PList_next(list, &temp)) {
      numPt += edge->getNbHighOrderPoints();
    }
    PList_delete(list);

    return numPt;
  }

  // -------------------------------------------------------------------
  //! Aligns the face with a set of points
  //! will return zero if not successful
  //! otherwise the absolute value will contain the index of the element vertex coinciding with the first vertex
  //! the sign indicates the orientation 
  int F_align(pFace pf,pVertex pv1,pVertex pv2,pVertex pv3,pVertex pv4)
  {
    return pf->align(pv1,pv2,pv3,pv4);
  }

  // -------------------------------------------------------------------
  //! Gets coordinates of vertices of the face, not including high-order points. \ingroup face
  void F_coordP1(const pFace face, double xyz[][3])
  {
    int iNode = 0;

    // Summits of the face (first order nodes)
    pPList fVerts = F_vertices(face,1);
    void * temp = NULL;
    while ( pVertex pV = (pVertex)PList_next(fVerts,&temp)) {
      V_coord(pV,xyz[iNode]);
      iNode++;
    }
    PList_delete(fVerts);
  }

  // -------------------------------------------------------------------
  //! Gets coordinates of points of the face including high-order points. \ingroup face
  void F_coord(const pFace face, double xyz[][3])
  {
    // Summits of the face (first order nodes)
    F_coordP1(face,xyz);
    int iNode = face->getNbNodes();

    // points on edges (higher order nodes)
    pPList fEdges = F_edges(face);
    void * temp = NULL;
    while ( pEdge edge = (pEdge) PList_next(fEdges, &temp)) {
      int nEPts = E_numPoints(edge);
      for (int iEP=0; iEP<nEPts; iEP++) {
        pPoint pP = E_point(edge,iEP);
        xyz[iNode][0] = P_x(pP); xyz[iNode][1] = P_y(pP); xyz[iNode][2] = P_z(pP);
        iNode++;
      }
    }
    PList_delete(fEdges);

    // points on face (higher order nodes)
    int nFPts = F_numPoints(face);
    for (int iFP=0; iFP<nFPts; iFP++) {
      pPoint pP = F_point(face,iFP);
      xyz[iNode][0] = P_x(pP); xyz[iNode][1] = P_y(pP); xyz[iNode][2] = P_z(pP);
      iNode++;
    }
  }

  // -------------------------------------------------------------------
  //! Gets the parametric coordinates of the summits of the face \ingroup face
  //! in the geometric entity on which it is classified.
  //! Returns false if parametric coordinates are not available.
  bool F_params(const pFace face, double u[][2])
  {
#ifdef _HAVE_GMSH_
    pGEntity faceGE = F_whatIn(face);
    int gDim = GEN_type(faceGE);
    if ( gDim != 2 ) return false;

    for (int iV=0; iV<F_numVertices(face); iV++)
      {
        pVertex pv = F_vertex(face,iV);
        pGEntity vG = EN_whatIn(pv);
        int vGDim = GEN_type(vG);
        
        switch ( vGDim ) {
        case 3: throw;
        case 2: {
          assert ( vG == faceGE );
          if ( !V_params(pv,&u[iV][0],&u[iV][1]) ) return false;
          break;
        }
        case 1: {
          double tmp0,tmp1;
          if ( !V_params(pv,&tmp0,&tmp1) ) return false;
          GE_reparamOnFace( (pGEdge)vG, (pGFace)faceGE, tmp0, u[iV], NULL );
          break;
        }
        case 0: {
          GV_reparamOnFace( (pGVertex)vG, (pGFace)faceGE, u[iV], NULL );
          break;
        }
        }
      }
    return true;
#else
    return false;
#endif
  }

  // -------------------------------------------------------------------
  //! Returns area of a triangular face \ingroup face
  //! if 'dir' is not NULL, the area is signed
  double F_area(const pFace face, const double * dir)
  {
    double fxyz[4][3];
    F_coordP1(face,fxyz);
    return XYZ_F_area(fxyz,dir);
  }

  // -------------------------------------------------------------------
  //! Returns area of a triangular face \ingroup face
  //! if 'dir' is not NULL, the area is signed
  double XYZ_F_area(const double xyz[][3], const double * dir)
  {
    double areaSq = XYZ_F_areaSq(xyz, dir);
    if ( areaSq >= 0. ) return sqrt(areaSq);
    return -1.0 * sqrt ( -1.0*areaSq );
  }

  // -------------------------------------------------------------------
  //! Returns square area of a triangular face. \ingroup face
  //! if 'dir' is not NULL, the square area is signed
  double F_areaSq(const pFace face, const double * dir)
  {
    double xyz[3][3];
    F_coordP1(face,xyz);
    return XYZ_F_areaSq(xyz,dir);
  }

  // -------------------------------------------------------------------
  //! Returns square area of a triangular face with coordinates xyz. \ingroup face
  //! if 'dir' is not NULL, the square area is signed
  double XYZ_F_areaSq(const double xyz[][3], const double * dir)
  {
    double e01[3], e02[3], normal[3];
    diffVec(xyz[1],xyz[0],e01);
    diffVec(xyz[2],xyz[0],e02);
    crossProd(e01,e02,normal);
    if( dir && ( dotProd(dir,normal) < 0. ) ) {
      return -0.25 * dotProd(normal,normal);
    }
    return 0.25 * dotProd(normal,normal);
  }
  // -------------------------------------------------------------------
  //! Returns the center of a triangular linear face. \ingroup face
  void F_center(const pFace face, double center[])
  {
    double fxyz[3][3];
    F_coordP1(face,fxyz);
    center[0] = MAdTHIRD * ( fxyz[0][0] + fxyz[1][0] + fxyz[2][0] );
    center[1] = MAdTHIRD * ( fxyz[0][1] + fxyz[1][1] + fxyz[2][1] );
    center[2] = MAdTHIRD * ( fxyz[0][2] + fxyz[1][2] + fxyz[2][2] );
  }
  
  // -------------------------------------------------------------------
  //! Returns the coordinates (u,v) of the point in the parent element
  void F_linearParams(const pFace pF, const double xyz[3], 
                      double res[2])
  {
    // using barycentric coordinates
    double fxyz[3][3];
    F_coordP1(pF,fxyz);
  
    Tri_linearParams(fxyz,xyz,res);
  }

  // -------------------------------------------------------------------
  //! Computes the normal vector of the face \ingroup face
  //! The normal is not normalized
  void F_normal(const pFace face, double normal[3])
  {
    double xyz[4][3];
    F_coordP1(face,xyz);
    double v01[3],v02[3];
    diffVec(xyz[1],xyz[0],v01);
    diffVec(xyz[2],xyz[0],v02);
    crossProd(v01,v02,normal);
  }

  // -------------------------------------------------------------------
  //! Computes the normal vector of a triangle with given coordinates \ingroup face
  //! The normal is not normalized
  void XYZ_F_normal(const double xyz[3][3], double normal[3])
  {
    double v01[3],v02[3];
    diffVec(xyz[1],xyz[0],v01);
    diffVec(xyz[2],xyz[0],v02);
    crossProd(v01,v02,normal);
  }
  
  // -------------------------------------------------------------------
  //! Computes the ratio between the regions around the face \ingroup face
  //! Returns false if there is less than 2 regions
  bool F_volumeRatio(const pFace face, double * ratio) 
  {
    const pRegion pr0 = F_region(face,0);
    const pRegion pr1 = F_region(face,1);
  
    if ( !( pr0 && pr1 ) ) return false;
  
    double vol0 = R_volume(pr0);
    double vol1 = R_volume(pr1);
  
    if( vol0 <= 0. || vol1 <= 0. ) *ratio = -1.;
    else *ratio = ( vol0 < vol1 ) ? vol1/vol0 : vol0/vol1;
  
    return true;
  }

  // -------------------------------------------------------------------
  //! Computes the maximal volume ratio among the faces \ingroup face
  double F_worstVolumeRatio(const std::set<pFace> faces)
  {
    double maxRatio = 0.;
  
    std::set<pFace>::const_iterator fIter = faces.begin();
    for(; fIter != faces.end(); fIter++ ) {
      double ratio;
      if ( !F_volumeRatio(*fIter,&ratio) ) continue;
      if( ratio < 0. )    return -1;
      if( ratio > maxRatio ) maxRatio = ratio;
    }
  
    return maxRatio;
  }

  // -------------------------------------------------------------------
  //! Returns the vertex of face opposite to edge \ingroup face
  pVertex F_edOpVt(const pFace face, const pEdge edge)
  {
    assert(face->getNbEdges() == 3);

    pVertex v0 = E_vertex(edge,0);
    pVertex v1 = E_vertex(edge,1);

    for (int iE=0; iE<3; iE++) {
      pEdge pE = face->getEdge(iE);
      if (pE==edge) continue;
      pVertex pEv0 = E_vertex(pE,0);
      if ( pEv0 != v0 && pEv0 != v1 ) return pEv0;
      else return E_vertex(pE,1);
    }

    return NULL;
  }

  // -------------------------------------------------------------------
  //! Returns the edge of face opposite to vertex \ingroup face
  pEdge F_vtOpEd(const pFace face, const pVertex vertex)
  {
    assert(face->getNbEdges() == 3);

    for (int iE=0; iE<3; iE++) {
      pEdge pE = face->getEdge(iE);
      if ( E_vertex(pE,0) != vertex &&
           E_vertex(pE,1) != vertex ) {
        return pE;
      }
    }
    return NULL;
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Edge operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  //! Dumps informations about the edge. \ingroup edge
  void E_info(const pEdge edge, std::string name, std::ostream& out)
  {
    out << "\n";
    out << "Edge \'" << name << "\' " << edge 
        << ", id: " << EN_id((pEntity)edge) << "\n";
    pGEntity pGE = EN_whatIn((pEntity)edge);
    out << "  - Classification: dim: " << GEN_type(pGE) << ", tag: " << GEN_tag(pGE) << "\n";
    pVertex pv0 = E_vertex(edge,0);
    pVertex pv1 = E_vertex(edge,1);
    out << "  - Vertices:\n"
        << "     * Vertex 0:\n";
    V_info(pv0, "", out);
    out << "     * Vertex 1:\n";
    V_info(pv1, "", out);
    out << "\n";
  }

  // -------------------------------------------------------------------
  //! Returns the number of regions attached to the edge pe \ingroup edge
  int E_numRegions(pEdge pe)
  {
    pPList rlist = E_regions(pe);
    int num = PList_size(rlist);
    PList_delete(rlist);
    return num;
  }

  // -------------------------------------------------------------------
  //! Returns the number of vertices in edge pe \ingroup edge
  //! Guess what, there are two !
  int E_numVertices(pEdge pe)
  {
    return 2;
  }

  // -------------------------------------------------------------------
  //! Returns the n-th vertex of edge pe \ingroup edge
  /*! \warning will throw an error for n>1 */
  pVertex E_vertex(pEdge pe, int n)
  {
    switch(n){
    case 0 : return pe->p1;
    case 1 : return pe->p2;
    }
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns a list of the vertices of edge pe \ingroup edge
  pPList E_vertices(const pEdge pe)
  {
    pPList vList;
    PList_append(vList,pe->p1);
    PList_append(vList,pe->p2);
    return vList;
  }

  // -------------------------------------------------------------------
  //! Returns the number of faces attached to edge pe \ingroup edge
  int E_numFaces(pEdge pe)
  {
    return pe->numfaces() ;
  }  

  // -------------------------------------------------------------------
  //! Returns a list of regions attached to edge pe \ingroup edge
  pPList E_regions(pEdge e)
  {
    pPList pl = PList_new();
    for (int i=0;i<e->numfaces();++i)
      {
        pFace pf = e->faces(i);
        for(int k=0;k<pf->getNbRegions();k++)
          PList_appUnique(pl,pf->getRegion(k));
      }
    return pl;
  }

  // -------------------------------------------------------------------
  //! Returns the number of higher-order points on edge pe, excluding vertices \ingroup edge
  int E_numPoints(pEdge pE)
  {
    return pE->getNbHighOrderPoints ();
  }

  // -------------------------------------------------------------------
  //! Returns the n-th face attached to edge pe \ingroup edge
  pFace E_face(pEdge pe, int n)
  {
    return pe->faces(n);
  }

  // -------------------------------------------------------------------
  //! Returns the list of faces connected to the edge \ingroup edge
  pPList E_faces(const pEdge pe)
  {
    pPList eFaces = PList_new();
    for (int iF=0; iF < pe->numfaces(); iF++) {
      PList_append(eFaces,pe->faces(iF));
    }
    return eFaces;    
  }

  // -------------------------------------------------------------------
  //! Returns 1 of ent is one of the vertices of edge pe, else 0 \ingroup edge
  int E_inClosure(pEdge pe, pEntity ent)
  {
    if (pe->p1 == ent) return 1;
    if (pe->p2 == ent) return 1;
    return 0;
  }

  // -------------------------------------------------------------------
  //! Returns the first other face attached to edge pe found \ingroup edge
  //! Warns if more than 2 faces around edge, returns NULL if no other face.
  pFace E_otherFace(pEdge edge, pFace face)
  {
    int nF = edge->numfaces();
    if ( nF > 2 ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                    "There are more than 2 faces around this edge, returning the first candidate");
    }

    for (int i=0;i<nF;i++)
      {
        pFace pf = edge->faces(i);
        if ( pf != face ) return pf;
      }
   
    return NULL;
  }

  // -------------------------------------------------------------------
  //! Returns the other face attached to edge pe belonging to the closure of region r \ingroup edge
  pFace E_otherFace(pEdge edge, pFace face, pRegion r)
  {
    for (int i=0;i<r->getNbFace();i++)
      {
        pFace f1 = r->getFace(i);
        if (f1 != face && F_inClosure (f1,edge) ) return f1;
      }
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns the other vertex of edge pe \ingroup edge
  pVertex  E_otherVertex(pEdge e, pVertex v)
  {
    if(e->p1 == v)return e->p2;
    if(e->p2 == v)return e->p1;
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns the n-th higher order point - excluding vertices - of edge pe \ingroup edge
  pPoint E_point(pEdge e, int n)
  {
    return e->getHighOrderPoint (n);
  }

  // -------------------------------------------------------------------
  //! Returns the geometrical entity on which entity 'pe' is classified \ingroup entity
  pGEntity EN_whatIn(pEntity pe)
  {
    return pe->g;
  }

  // -------------------------------------------------------------------
  //! Returns the type/dimension of the geometrical entity on which edge pe is classified \ingroup edge
  pGEntity  E_whatIn(pEdge pe)
  {
    return pe->g;
  }

  // -------------------------------------------------------------------
  //! Classify entity e on the geometrical entity 'what'  \ingroup entity
  void EN_setWhatIn(pEntity e, pGEntity what){e->g = what;}

  // -------------------------------------------------------------------
  //! Classify edge e on the geometrical entity 'what'  \ingroup entity
  void E_setWhatIn(pEdge edge, pGEntity what)
  {
    edge->g = what;
    pPoint p1 = edge->p1;
    if ( GEN_type(p1->g) > GEN_type(what) ) V_setWhatIn(p1,what);
    pPoint p2 = edge->p2;
    if ( GEN_type(p2->g) > GEN_type(what) ) V_setWhatIn(p2,what);
  }

  // -------------------------------------------------------------------
  //! Returns the type/dimension of the geometrical entity on which edge pe is classified \ingroup edge
  int E_whatInType(pEdge e)   { return EN_whatInType(e);}

  // -------------------------------------------------------------------
  //! Aligns the current edge with the node order \ingroup edge
  //! returns 1 if the edge is already aligned, -1 if the edge is inverted and 0 if the edge does not correspond
  int E_align(pEdge pe,pVertex pv1,pVertex pv2) {return pe->align(pv1,pv2);}

  // -------------------------------------------------------------------
  //! Returns 1 of the edge is oriented from first vertex \ingroup edge
  //! to second one, 0 otherwise
  int E_dir(pEdge pe, pVertex pv1, pVertex pv2) {
    if (pe->p1 == pv1 && pe->p2 == pv2) return 1;
    if (pe->p1 == pv2 && pe->p2 == pv1) return 0;
    throw;
  }

  // -------------------------------------------------------------------
  //! Returns coordinates of vertices of the edge not including high-order points. \ingroup edge
  void E_coordP1(const pEdge edge, double xyz[][3]) 
  {
    // points of extremities
    MDB_Point * pP1 = edge->p1;
    xyz[0][0] = P_x(pP1);
    xyz[0][1] = P_y(pP1);
    xyz[0][2] = P_z(pP1);
    MDB_Point * pP2 = edge->p2;
    xyz[1][0] = P_x(pP2);
    xyz[1][1] = P_y(pP2);
    xyz[1][2] = P_z(pP2);
  }

  // -------------------------------------------------------------------
  //! Returns coordinates of vertices of the edge including high-order points. \ingroup edge
  void E_coord(const pEdge edge, double xyz[][3]) 
  {
    E_coordP1(edge,xyz);

    // points on edges (higher order nodes)
    int n = E_numPoints(edge);
    for (int i=0; i<n; i++) {
      pPoint pnt = E_point(edge,i);
      xyz[2+i][0] = P_x(pnt); xyz[2+i][1] = P_y(pnt); xyz[2+i][2] = P_z(pnt);
    }
  }

  // -------------------------------------------------------------------
  //! Gets the parametric coordinates of the summits of the edge \ingroup edge
  //! in the geometric entity on which it is classified.
  //! Returns false if parametric coordinates are not available.
  bool E_params(const pEdge edge, double u[2][2])
  {
#ifdef _HAVE_GMSH_

    pGEntity edgeGE = E_whatIn(edge);
    int gDim = GEN_type(edgeGE);
  
    if ( gDim != 1 && gDim != 2 ) return false;

    for (int iV=0; iV<2; iV++) {
    
      pVertex pv = E_vertex(edge,iV);
      pGEntity vG = EN_whatIn(pv);
      int vGDim = GEN_type(vG);

      // --------------------------------------------------------------------
      // edge classified on a surface: we want the parameters on the surface
      // --------------------------------------------------------------------
      if ( gDim == 2 ) {

        switch ( vGDim ) {
        case 3: throw;
        case 2: {
          if ( !V_params(pv,&u[iV][0],&u[iV][1]) ) return false;
//           printf("node on face %d: %lf %lf\n",GEN_tag(edgeGE),u[iV][0],u[iV][1]);
          break;
        }
        case 1: {
          double tmp0,tmp1;
          if ( !V_params(pv,&tmp0,&tmp1) ) return false;

          double * otherU = NULL;
          if ( GE_isSeam( (pGEdge)vG, (pGFace)edgeGE ) )
            {
              pVertex otherV = E_vertex(edge,1-iV);
              double tmp20,tmp21;
              if ( V_whatInType(otherV)==2 && V_params(otherV,&tmp20,&tmp21) ) {
                otherU = new double[2];
                otherU[0] = tmp20;
                otherU[1] = tmp21;
              }
              else if ( V_whatInType(otherV)==1 && V_params(otherV,&tmp20,&tmp21) ) {
                if ( GE_isSeam( (pGEdge)V_whatIn(otherV), (pGFace)edgeGE ) ) {
                  MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                                "Found a surface with 2 seams");
                  return false;
                }
                otherU = new double[2];
                GE_reparamOnFace( (pGEdge)V_whatIn(otherV), (pGFace)edgeGE, 
                                  tmp20, &(otherU[0]) );
                otherU[1] = -1.;
              }
            }
          GE_reparamOnFace( (pGEdge)vG, (pGFace)edgeGE, tmp0, u[iV], otherU );
//           printf("node from edge %d on face %d: %lf %lf\n",GEN_tag(vG),GEN_tag(edgeGE),u[iV][0],u[iV][1]);
          break;
        }
        case 0: {
          double tmp0,tmp1;
          if ( !V_params(pv,&tmp0,&tmp1) ) return false;

          double * otherU = NULL;
          if ( GV_isOnSeam( (pGVertex)vG, (pGFace)edgeGE ) ) 
            {
              pVertex otherV = E_vertex(edge,1-iV);
              double tmp20,tmp21;
              if ( V_whatInType(otherV)==2 && V_params(otherV,&tmp20,&tmp21) ) {
                otherU = new double[2];
                otherU[0] = tmp20;
                otherU[1] = tmp21;
              }
              else if ( V_whatInType(otherV)==1 && V_params(otherV,&tmp20,&tmp21) ) {
                if ( GE_isSeam( (pGEdge)V_whatIn(otherV), (pGFace)edgeGE ) ) {
                  MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                                "Found a surface with 2 seams");
                  return false;
                }
                otherU = new double[2];
                GE_reparamOnFace( (pGEdge)V_whatIn(otherV), (pGFace)edgeGE, 
                                  tmp20, &(otherU[0]) );
                otherU[1] = -1.;
              }
            }
          GV_reparamOnFace( (pGVertex)vG, (pGFace)edgeGE, u[iV], otherU );
//           printf("vert %d on face %d: %lf %lf\n",GEN_tag(vG),GEN_tag(edgeGE),u[iV][0],u[iV][1]);
          break;
        }
        }
      }

      // --------------------------------------------------------------
      // edge classified on a line: we want the parameters on the line
      // --------------------------------------------------------------
      else if ( gDim == 1) {

        switch ( vGDim ) {
        case 3: throw;
        case 2: throw;
        case 1: {
          double tmp;
          if ( !V_params(pv,&u[iV][0],&tmp) ) return false;
          break;
        }
        case 0: {
          if ( !(pv->isParametric() ) ) return false;

          double otherU = -1.;
          pVertex otherV = E_vertex(edge,1-iV);
          double tmp20,tmp21;
          if ( V_whatInType(otherV)==1 && V_params(otherV,&tmp20,&tmp21) ) {
            otherU = tmp20;
          }
          else if ( V_whatInType(otherV)==0 ) {
            GV_reparamOnEdge( (pGVertex)V_whatIn(otherV), (pGEdge)edgeGE, &otherU );
          }

          GV_reparamOnEdge( (pGVertex)vG, (pGEdge)edgeGE, &u[iV][0], otherU );
          break;
        }
        }

      }

      else throw;
    }

    return true;
#else
    return false;
#endif
  }

  // -------------------------------------------------------------------
  //! Returns the physical length of the edge \ingroup edge
  double E_length(const pEdge edge)
  {
    double eCoords[2][3];
    E_coordP1(edge,eCoords);
    double e[3];
    diffVec(eCoords[1],eCoords[0],e);
    return sqrt ( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
  }

  // -------------------------------------------------------------------
  //! Returns the physical square length of the edge \ingroup edge
  double E_lengthSq(const pEdge edge)
  {
    double eCoords[2][3];
    E_coordP1(edge,eCoords);
    double e[3];
    diffVec(eCoords[1],eCoords[0],e);
    return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
  }
  
  // -------------------------------------------------------------------
  double E_linearParams(const pEdge pE, 
                        const pVertex pv)
  {
    double xyz[3];
    V_coord(pv,xyz);
    return E_linearParams(pE,xyz);
  }
  
  // -------------------------------------------------------------------
  double E_linearParams(const pEdge pE, 
                        const double xyz[3])
  {
    double eCoords[2][3];
    E_coordP1(pE,eCoords);

    double v01[3];
    diffVec(eCoords[1],eCoords[0],v01);

    double vX1[3];
    diffVec(eCoords[1],xyz,vX1);
    double temp = dotProd(v01,vX1);
    if ( temp < MAdTOL ) return 1.;

    double v0X[3];
    diffVec(xyz,eCoords[0],v0X);
    double temp2 = dotProd(v01,v0X);
    if ( temp2 < MAdTOL ) return 0.;

    double ratio = temp2 / temp;
    return ( ratio / (1. + ratio ) );
  }

  // -------------------------------------------------------------------
  //! Returns the center of a linear edge. \ingroup edge
  void E_center(const pEdge edge, double center[3])
  {
    double exyz[2][3];
    E_coordP1(edge,exyz);
    center[0] = 0.5 * ( exyz[0][0] + exyz[1][0] );
    center[1] = 0.5 * ( exyz[0][1] + exyz[1][1] );
    center[2] = 0.5 * ( exyz[0][2] + exyz[1][2] );
  }

  // -------------------------------------------------------------------
  //! Returns the center of the cavity surrounding an edge. \ingroup edge
  void E_cavityCenter(const pEdge edge, double center[3])
  {
    if( E_whatInType(edge)!=3 ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "not implemented for edges classified on <3D");
    }

    center[0] = 0.; center[1] = 0.; center[2] = 0.;

    int nbF = 0;
    pPList vFaces = E_faces(edge);
    void * temp = NULL;
    while ( pFace face = (pFace)PList_next(vFaces,&temp) ) {
      pVertex oppV = F_edOpVt(face,edge);
      double xyz[3];
      V_coord(oppV,xyz);
      center[0] += xyz[0]; center[1] += xyz[1]; center[2] += xyz[2];
      nbF++;
    }
    PList_delete(vFaces);

    double invNbF = 1. / nbF;
    center[0] *= invNbF; center[1] *= invNbF; center[2] *= invNbF;
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Vertex operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  //! Dumps informations about the vertex. \ingroup vertex
  void V_info(const pVertex vertex, std::string name, std::ostream& out)
  {
    out << "\n";
    out << "Vertex \'" << name << "\' " << vertex
        << ", id: " << EN_id((pEntity)vertex) << "\n";
    pGEntity pGE = EN_whatIn((pEntity)vertex);
    out << "  - Classification: dim: " << GEN_type(pGE) << ", tag: " << GEN_tag(pGE) << "\n";
    out << "  - Num edges: " << vertex->edges.size() << "\n";
    out << "  - Parameters:\n";
    double u,v;
    if ( vertex->getParams(&u,&v) ) {
      out << "     Point is parametric: (u,v) = ("<<u<<", "<<v<<")\n";
    }
    else {
      out <<"      Point is not parametric\n";
    }
    out << "  - Coordinates:\n";
    double xyz[3];
    V_coord(vertex,xyz);
    out << "     ( "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<" )\n";
    out << "\n";
    
  }

  // -------------------------------------------------------------------
  //! Returns the type/dimension of the geometrical entity on which vertex e is classified \ingroup vertex
  int V_whatInType(pVertex e)
  {
    return EN_whatInType(e);
  }

  // -------------------------------------------------------------------
  //! Classify vertex e on the geometrical entity 'what'  \ingroup vertex
  void V_setWhatIn(pVertex e, pGEntity what)
  {
    e->g = what;
  }

  // -------------------------------------------------------------------
  //! Returns the geometrical entity on which vertex pv is classified \ingroup vertex
  pGEntity V_whatIn(pVertex pv)
  { 
    return pv->g;
  }

  // -------------------------------------------------------------------
  //! Returns the list of regions attached to vertex v \ingroup vertex
  pPList V_regions(pVertex v)
  {
    //  3D SPECIFIC printf("coucou\n");
    pPList pl = PList_new();
    for (unsigned int i=0;i<v->edges.size();++i)
      {
        pEdge pe = v->edges[i];
        for (int j=0;j<pe->numfaces();j++)
          {
            pFace pf = pe->faces(j);
            for(int k=0;k<pf->getNbRegions();k++)
              PList_appUnique(pl,pf->getRegion(k));
          }
      }  
    return pl;
  }

  // -------------------------------------------------------------------
  //! Copies the coordinates of vertex pv in point \ingroup vertex
  pPoint  V_point(pVertex pv)
  {
    return pv;
  }

  // -------------------------------------------------------------------
  //! Returns number of edges attached to vertex pv \ingroup vertex
  int V_numEdges(pVertex  pv)
  {
    return pv->edges.size();
  }

  // -------------------------------------------------------------------
  //! Returns number of faces attached to vertex pv \ingroup vertex
  int V_numFaces(pVertex  pv)
  {
    pPList faces =  V_faces(pv);
    int num = PList_size(faces);
    PList_delete(faces);
    return num;
  }

  // -------------------------------------------------------------------
  //! Returns number of regions attached to vertex pv \ingroup vertex
  int V_numRegions(pVertex  pv)
  {
    pPList regions =  V_regions(pv);
    int num = PList_size(regions);
    PList_delete(regions);
    return num;
  }

  // -------------------------------------------------------------------
  //! Returns n-th edge attached to vertex pv \ingroup vertex
  pEdge V_edge(pVertex pv, int n)
  {
    return pv->edges[n];
  }


  // -------------------------------------------------------------------
  //! Returns the identity tag of the vertex pv \ingroup vertex
  int V_id(const pVertex pv) {
    return pv->iD;
  }

  // -------------------------------------------------------------------
  //! Copy position of vertex p to xyz \ingroup vertex
  void V_coord(const pVertex p, double xyz[3])
  {
    xyz[0] = p->X;
    xyz[1] = p->Y;
    xyz[2] = p->Z;
  }

  // -------------------------------------------------------------------
  //! Gets the parametric coordinates of the vertex. \ingroup vertex
  //! Returns false is the vertex is not parametric. 
  bool V_params(pVertex pv, double * u, double * v)
  {
    return pv->getParams(u,v);
  }

  // -------------------------------------------------------------------
  //! Return the list of faces attached to vertex v \ingroup vertex
  pPList V_faces(pVertex v)
  {
    pPList pl = PList_new();
    for (unsigned int i=0;i<v->edges.size();++i)
      {
        pEdge pe = v->edges[i];
        for (int j=0;j<pe->numfaces();j++)
          {
            pFace pf = pe->faces(j);
            PList_appUnique(pl,pf);
          }
      }
  
    return pl;
  }

  // -------------------------------------------------------------------
  //! Return the list of edges attached to vertex v \ingroup vertex
  pPList V_edges(pVertex v)
  {
    pPList pl = PList_new();
    for (unsigned int i=0; i<v->edges.size(); ++i) {
      pEdge pe = v->edges[i];
      PList_append(pl,pe);
    }
  
    return pl;
  }

  // -------------------------------------------------------------------
  //! Compute the mean length square of the edges adjacent to the vertex \ingroup vertex
  double V_meanEdgeLenSq(const pVertex pv)
  {
    double lenSq = 0.;
    for( int iE=0; iE<V_numEdges(pv); iE++ ) {
      lenSq += E_lengthSq( V_edge(pv,iE) );
    }
    lenSq /= V_numEdges(pv);
    return lenSq;
  }

  // -------------------------------------------------------------------
  //! Returns the center of the cavity surrounding a vertex. \ingroup vertex
  void V_cavityCenter(const pVertex vertex, double center[3])
  {
    center[0] = 0.; center[1] = 0.; center[2] = 0.;

    int nbE = 0;
    pPList vEdges = V_edges(vertex);
    void * temp = NULL;
    while ( pEdge edge = (pEdge)PList_next(vEdges,&temp) ) {
      pVertex otherV = E_otherVertex(edge,vertex);
      double xyz[3];
      V_coord(otherV,xyz);
      center[0] += xyz[0]; center[1] += xyz[1]; center[2] += xyz[2];
      nbE++;
    }
    PList_delete(vEdges);

    double invNbE = 1. / nbE;
    center[0] *= invNbE; center[1] *= invNbE; center[2] *= invNbE;
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Point operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  //! Returns x coordinate of point p \ingroup point
  double P_x(pPoint p){return p->X;}

  // -------------------------------------------------------------------
  //! Returns y coordinate of point p \ingroup point
  double P_y(pPoint p){return p->Y;}

  // -------------------------------------------------------------------
  //! Returns z coordinate of point p \ingroup point
  double P_z(pPoint p){return p->Z;}

  // -------------------------------------------------------------------
  //! Deletes point p \ingroup point
  void P_delete(pPoint p)
  {
    delete p;
  }

  // -------------------------------------------------------------------
  //! Sets the coordinates of point p \ingroup point
  void P_setPos(pPoint p , double x, double y, double z)
  {
    p->X = x;
    p->Y = y;
    p->Z = z;
  }

  // -------------------------------------------------------------------
  //! Sets the identity of point p \ingroup point
  void P_setID(pPoint p, int id){p->iD = id;}

  // -------------------------------------------------------------------
  //! Returns the identity of point p \ingroup point
  int P_id(pPoint p){return p->iD;}

  // -------------------------------------------------------------------
  double P_param1(pPoint p)
  {
    pMeshDataId id = MD_lookupMeshDataId("_param1");
    double pp;
    EN_getDataDbl(p, id,&pp);
    return pp;
  }

  // -------------------------------------------------------------------
  void P_setParam1(pPoint p, double param)
  {
    pMeshDataId id = MD_lookupMeshDataId("_param1");
    EN_attachDataDbl(p, id, param);
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Entity operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  //! Print informations about the entity
  void EN_info(const pEntity e, std::string name, std::ostream& out)
  {
    switch ( EN_type(e) ) {
    case 3: R_info((pRegion)e,name,out);
    case 2: F_info((pFace)e,  name,out);
    case 1: E_info((pEdge)e,  name,out);
    case 0: V_info((pVertex)e,name,out);
    }
  }

  // -------------------------------------------------------------------
  //! Returns the identity tag for entity e \ingroup entity
  int EN_id(pEntity pe)
  {
    return pe->iD;
  }

  // -------------------------------------------------------------------
  //! Sets the identity tag for entity e \ingroup entity
  void EN_setID(pEntity pe, int id){
    pe->iD=id;
  }

  // -------------------------------------------------------------------
  //! Returns the type(dimension) of the geometrical entity on which entity e is classified \ingroup entity
  int EN_whatInType(pEntity e){return GEN_type(e->g);}

  // -------------------------------------------------------------------
  //! Returns the number of vertices in the entity \ingroup entity
  int EN_numVertices(const pEntity e)
  {
    switch ( EN_type(e) ) {
    case 3: return R_numVertices((pRegion)e);
    case 2: return F_numVertices((pFace)e);
    case 1: return E_numVertices((pEdge)e);
    case 0: return 1;
    }
    return -1;
  }

  // -------------------------------------------------------------------
  //! Returns the list of the vertices of the entity \ingroup entity
  pPList EN_vertices(const pEntity e)
  {
    pPList vList;
    switch ( EN_type(e) ) {
    case 3: { vList = R_vertices((pRegion)e); break; }
    case 2: { vList = F_vertices((pFace)e,1); break; }
    case 1: { vList = E_vertices((pEdge)e); break; }
    case 0: { vList = PList_new(); PList_append(vList,e); break; }
    }
    return vList;
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  /* mesh entity creation routines */

  //! Creates a region, classified on gent, in mesh m using specified faces \ingroup mesh
  pRegion M_createR(pMesh mesh, int nFace, pFace faces[], pGEntity gent)
  {
    pRegion pr(0);
    switch(nFace)
      {
      case 4 : pr =  mesh->add_tet((MDB_Triangle*)faces[0], (MDB_Triangle*)faces[1],
                                   (MDB_Triangle*)faces[2], (MDB_Triangle*)faces[3],
                                   gent);
        break;
      case 5 : pr =  mesh->add_prism((MDB_Triangle*)faces[0], (MDB_Triangle*)faces[1],
                                     (MDB_Quad*)faces[2], (MDB_Quad*)faces[3],
                                     (MDB_Quad*)faces[4], gent);
        break;
      case 6 : pr =  mesh->add_hex((MDB_Quad*)faces[0], (MDB_Quad*)faces[1],
                                   (MDB_Quad*)faces[2], (MDB_Quad*)faces[3],
                                   (MDB_Quad*)faces[4], (MDB_Quad*)faces[5],gent);
        break;      
      }
    return pr;
  }

  //! Creates a region, classified on gent, in mesh m using specified vertices \ingroup mesh
  pRegion M_createR(pMesh mesh, int nVert, pVertex verts[], pGEntity gent)
  {
    pRegion pr(0);
    switch(nVert)
      {
      case 4 : pr =  mesh->add_tet(verts[0], verts[1], verts[2], verts[3],
                                   gent);
        break;
      case 6 : pr =  mesh->add_prism(verts[0], verts[1], verts[2], 
                                     verts[3], verts[4], verts[5],
                                     gent);
        break;
      case 8 : pr =  mesh->add_hex(verts[0], verts[1], verts[2], verts[3],
                                   verts[4], verts[5], verts[6], verts[7],
                                   gent);
        break;      
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Not a valid number of vertices: %d",
                                    nVert);
      }
    return pr;
  }

  //! Creates a region, classified on gent, in mesh m using specified vertices \ingroup mesh
  pRegion M_createR(pMesh mesh, int nVert, int vId[], pGEntity gent)
  {
    pRegion pr(0);
    switch(nVert)
      {
      case 4 : pr =  mesh->add_tet(vId[0], vId[1], vId[2], vId[3],
                                   gent);
        break;
      case 6 : pr =  mesh->add_prism(vId[0], vId[1], vId[2], 
                                     vId[3], vId[4], vId[5],
                                     gent);
        break;
      case 8 : pr =  mesh->add_hex(vId[0], vId[1], vId[2], vId[3],
                                   vId[4], vId[5], vId[6], vId[7],
                                   gent);
        break;      
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Not a valid number of vertices: %d",
                                    nVert);
      }
    return pr;
  }

  //! Creates a high order region, classified on gent, in mesh m using 
  //! specified vertices \ingroup mesh
  pRegion M_createTet(pMesh m, int order, bool serendip, int vId[], pGEntity gent)
  {
    return m->add_tet(gent,order,serendip,vId); 
  }

  //! Creates a face, classified on gent, in mesh m using specified edges \ingroup mesh
  /*! \warning currently only implemented for triangles */
  pFace M_createF(pMesh mesh, int nEdge, pEdge edges[], pGEntity gent)
  {
    pFace pf(0);
    switch(nEdge)
      {
      case 3: pf = mesh->add_triangle (edges[0], edges[1], edges[2], gent);
        break;
      case 4: pf = mesh->add_quad (edges[0], edges[1], edges[2], edges[3], gent);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Not a valid number of edges: %d", nEdge);
      }
    return pf;
  }

  //! Creates a face, classified on gent, in mesh m using specified vertices \ingroup mesh
  pFace M_createF(pMesh mesh, int nVert, pVertex verts[], pGEntity gent)
  {
    pFace pf(0);
    switch(nVert)
      {
      case 3: pf = mesh->add_triangle (verts[0], verts[1], verts[2], gent);
        break;
      case 4: pf = mesh->add_quad (verts[0], verts[1], verts[2], verts[3], gent);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Not a valid number of vertices: %d", nVert);
      }
    return pf;
  }

  //! Creates a face, classified on gent, in mesh m using specified vertices \ingroup mesh
  pFace M_createF(pMesh mesh, int nVert, int vId[], pGEntity gent)
  {
    pFace pf(0);
    switch(nVert)
      {
      case 3: pf = mesh->add_triangle (vId[0], vId[1], vId[2], gent);
        break;
      case 4: pf = mesh->add_quad (vId[0], vId[1], vId[2], vId[3], gent);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Not a valid number of vertices: %d", nVert);
      }
    return pf;
  }

  //! Creates a complete high order triangle, classified on gent, in mesh m using 
  //! specified vertices \ingroup mesh
  //! template is the following v0 ho[v0->v1] v1 ho[v1->v2] v2 ho[v2->v0]
  pFace M_createTri(pMesh mesh, int order, int vId[], pGEntity gent)
  {
    pFace pf(0);
    switch(order)
      {
      case 1: pf = mesh->add_triangle (vId[0], vId[1], vId[2], gent);
        break;
      case 2: pf = mesh->add_triangle (order, true, gent,
                                       vId[0], vId[1], vId[2], 
                                       vId[3], vId[4], vId[5] );
        break;
      case 3: pf = mesh->add_triangle (order, true, gent,
                                       vId[0], vId[1], vId[2], 
                                       vId[3], vId[4], vId[5],
                                       vId[6], vId[7], vId[8], vId[9] );
        break;
      case 4: pf = mesh->add_triangle (order, true, gent,
                                       vId[0], vId[1], vId[2], 
                                       vId[3], vId[4], vId[5],
                                       vId[6], vId[7], vId[8], vId[9],
                                       vId[10], vId[11], vId[12], vId[13], vId[14] );
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Complete triangles not defined for order %d",
                                    order);
      }
    return pf;

  }

  //! Create a first order edge, classified on ent, in mesh 'mesh' using 
  //! specified vertices \ingroup mesh
  pEdge M_createE(pMesh mesh, pVertex v1, pVertex v2, pGEntity ent)
  {
    return mesh->add_edge (v1,v2,ent);
  }

  //! Create a first order edge, classified on ent, in mesh 'mesh' using 
  //! specified vertices \ingroup mesh
  pEdge M_createE(pMesh mesh, int id1, int id2, pGEntity ent)
  {
    return mesh->add_edge (id1,id2,ent);
  }

  //! Create a high order edge, classified on ent, in mesh 'mesh' using 
  //! specified vertices \ingroup mesh
  pEdge M_createE(pMesh mesh, pVertex v1, pVertex v2, int order, 
                  pVertex internalPoints[], pGEntity ent)
  {
    return mesh->add_edge(v1,v2,ent,order,internalPoints);
  }

  pVertex M_createV(pMesh mesh, double x, double y, double z,
                    int patch, pGEntity ent)
  {
    return mesh->add_point (patch,x,y,z,ent);
  }
  pVertex M_createV2(pMesh mesh, double xyz[3],int patch, pGEntity ent)
  {
    return M_createV(mesh,xyz[0],xyz[1],xyz[2],patch,ent);
  }
  pVertex M_createVP(pMesh mesh, double x, double y, double z,
                     double u, double v, int patch, pGEntity ent)
  {
    return mesh->add_pointParam (patch,x,y,z,u,v,ent);
  }
  /* mesh entity deletion routines */

  //! Removes a region from the mesh, deletes it \ingroup mesh
  void M_removeRegion(pMesh m, pRegion region)
  {
    if (dynamic_cast<MDB_Tet*>(region))  m->del_tet((MDB_Tet*)region);
    else if (dynamic_cast<MDB_Hex*>(region))  m->del_hex((MDB_Hex*)region); 
    else if (dynamic_cast<MDB_Prism*>(region))  m->del_prism((MDB_Prism*)region); 
  }

  //! Removes a face from the mesh, deletes it \ingroup mesh
  /*! \warning currently only implemented for triangles */
  void M_removeFace(pMesh m, pFace face)
  {
    m->del_triangle((MDB_Triangle*)face);
  }

  //! Removes an edge from the mesh, deletes it \ingroup mesh
  void M_removeEdge(pMesh m, pEdge edge)
  {
    m->del_edge (edge);
  }

  //! Removes a vertex from the mesh, deletes it \ingroup mesh
  void M_removeVertex(pMesh m, pVertex vertex)
  {
    m->del_point(vertex);
  }

  //! Return edge between two vertices, if not found returns 0  \ingroup mesh
  pEdge E_exist(pVertex v1, pVertex v2)
  {
    MDB_VectorE edges = v1->edges;
    MDB_VectorE::iterator it  = edges.begin();
    MDB_VectorE::iterator ite = edges.end();
    while(it!=ite){
      pVertex p = (*it)->othervertex(v1);
      if(p==v2) return((*it));
      ++it;
    }
    return 0;
  }

  //! Returns triangle defined by an edge and a vertex. Returns 0 if not found.  \ingroup mesh
  pFace F_exist(pEdge edge, pVertex vertex)
  {
    for (int iF=0; iF<edge->numfaces(); iF++) {
      MDB_Face * face = edge->faces(iF);
      if(face->getNbEdges()!=3)
        continue;
      MDB_Triangle *tri=(MDB_Triangle*) face;
      MDB_Edge * otherE = tri->e1;
      if ( otherE == edge) otherE = tri->e2;
      if ( otherE->p1 == vertex || otherE->p2 == vertex ) {
        if ( edge->p1 == vertex || edge->p2 == vertex ) {
          F_info(tri);
          MAdMsgSgl::instance().error(__LINE__,__FILE__,"Bad triangle or bad input (edge: %p, vertex: %p)",
                                      edge,vertex);
        }
        return tri;
      }
    }
    return NULL;
  }

  //! Returns face defined by at most 4 entities, return 0 on failure  \ingroup mesh
  /*! entities e1..4 may be either vertices or edges */
  pFace F_exist(pEntity e1, pEntity e2, pEntity e3, pEntity e4){
    int typeEnt = EN_type(e1);
    switch(typeEnt) {
    case 0:
      { pVertex p1 = (pVertex) e1;
        pVertex p2 = (pVertex) e2;
        pVertex p3 = (pVertex) e3;
        pVertex p4 = (pVertex) e4;
        MDB_ListFace listFaces;
        p1->getFaces(listFaces);
        MDB_ListFace::iterator it  = listFaces.begin();
        MDB_ListFace::iterator ite = listFaces.end();
        while(it!=ite){
          pVertex p[4]={NULL,NULL,NULL,NULL};
          (*it)->getNodes(p);
          if (   (p[0]==p1 || p[0]==p2 || p[0]==p3 || p[0]==p4)
              && (p[1]==p1 || p[1]==p2 || p[1]==p3 || p[1]==p4)
              && (p[2]==p1 || p[2]==p2 || p[2]==p3 || p[2]==p4)
              && (p[3]==p1 || p[3]==p2 || p[3]==p3 || p[3]==p4))
                return (*it);
          ++it;
        }      
        return 0;
      }
      break;
    case 1:
      {
        pEdge ped1 = (pEdge) e1;
        pEdge ped2 = (pEdge) e2;
        pEdge ped3 = (pEdge) e3;
        pEdge ped4 = (pEdge) e4;
        for(int k=0 ; k<ped1->numfaces() ; k++) {
          pFace pface = ped1->faces(k);
          if ((pface->getEdge(0)==ped1 || pface->getEdge(1)==ped1 || pface->getEdge(2)==ped1  || pface->getEdge(3)==ped1)
            &&(pface->getEdge(0)==ped2 || pface->getEdge(1)==ped2 || pface->getEdge(2)==ped2  || pface->getEdge(3)==ped2)
            &&(pface->getEdge(0)==ped3 || pface->getEdge(1)==ped3 || pface->getEdge(2)==ped3  || pface->getEdge(3)==ped3)
            &&(pface->getEdge(0)==ped4 || pface->getEdge(1)==ped4 || pface->getEdge(2)==ped4  || pface->getEdge(3)==ped4))
              return (pface);
        }
        return 0;
      }
      break;
    default:
      throw;  
    }
    return 0;
  }

  /* deprecated (rather use ModelInterface) */
  void M_setGeomFeature(pMesh m, int tag, pGEntity geom)
  {
    bool found = false;
    std::multimap<int, pGEntity>::iterator it;
    for (it  = m->geomFeatures_Tags.lower_bound(tag);
         it != m->geomFeatures_Tags.upper_bound(tag);++it)
      if (it->second == geom) found = true;
    if (!found)
      m->geomFeatures_Tags.insert(std::pair<int,pGEntity>(tag, geom));
  }

  //! Returns the number of geometric entities associated to the mesh \ingroup mesh

  int    M_numGeomFeatures (pMesh pm)
  {
    return pm->geomFeatures_Tags.size(); 
  }

  //! Returns the id tag of the n-th geometric entity associated to the mesh \ingroup mesh
  int    M_geomFeatureId   (pMesh pm, int n)
  {
    int count = 0;
    for (std::multimap<int, pGEntity>::iterator it = pm->geomFeatures_Tags.begin();
         it != pm->geomFeatures_Tags.end();
         ++it)
      {
        if (count++ == n) return it->first;
      }
    throw;
  }

  //! Returns the number of geometric entities associated to the mesh with tag id \ingroup mesh
  /*! \warning id's are only guaranteed to be unique per type of geometric feature */
  pPGList M_geomFeature     (pMesh pm, int id)
  {
    pPGList pl = PGList_new();
    for (std::multimap<int, pGEntity>::iterator it = pm->geomFeatures_Tags.lower_bound(id) ;
         it != pm->geomFeatures_Tags.upper_bound(id) ; ++it)
      PGList_append(pl,it->second);
    return pl;  
  }

  //! Print a list of geometric entities attached to the mesh to standard output \ingroup mesh
  /*! \warning id's  are only guaranteed to be unique per type of geometric feature */
  void M_printGeomFeatures (pMesh pm)
  {
    std::cout << "\nPrinting geometric features: \n\n";
    std::multimap<int, pGEntity>::iterator it = pm->geomFeatures_Tags.begin();
    std::multimap<int, pGEntity>::iterator itEnd = pm->geomFeatures_Tags.end();
    for (;it != itEnd; it++) {
      std::cout << "geometric id: " << it->first << " geometric entity: " << it->second << "\n";
    }
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // PList operators
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  PList * PList_new()
  {
    return new PList();
  }

  // -------------------------------------------------------------------
  void PList_delete(PList * lst)
  {
    if (lst) delete lst;
    lst = NULL;
  }

  // -------------------------------------------------------------------
  void PList_clear(PList * lst)
  {
    lst->clear();
  }

  // -------------------------------------------------------------------
  PList * PList_appUnique(PList * lst, MDB_MeshEntity * ent)
  {
    for ( unsigned int i=0; i<lst->entities.size(); i++ ) {
      if ( lst->entities[i] == ent ) return lst;
    }
    lst->entities.push_back(ent);
    return lst;
  }

  // -------------------------------------------------------------------
  PList * PList_appPListUnique(PList * lst, PList * source)
  {
    for ( unsigned int iSrc=0; iSrc<source->entities.size(); iSrc++ ) {
      PList_appUnique(lst,source->entities[iSrc]);
    }
    return lst;
  }

  // -------------------------------------------------------------------
  PList * PList_append(PList * lst, MDB_MeshEntity * ent)
  {
    lst->entities.push_back(ent);
    return lst;
  }

  // -------------------------------------------------------------------
  int PList_size(PList * lst)
  {
    return lst->entities.size();
  }

  // -------------------------------------------------------------------
  MDB_MeshEntity * PList_item(PList * lst, int i)
  {
    return lst->entities[i];
  }

  // -------------------------------------------------------------------
  MDB_MeshEntity * PList_next(PList * lst, void **restart)
  {
    if( *(int*)(restart) >= (int)lst->entities.size() ) return NULL;
    return lst->entities[(*(int*)(restart))++];
  }

  // -------------------------------------------------------------------
  int PList_inList(PList * lst, MDB_MeshEntity * ent)
  {
    for ( unsigned int i=0; i<lst->entities.size(); i++ ) {
      if ( lst->entities[i] == ent ) return 1;
    }
    return 0;
  }

  // -------------------------------------------------------------------
  void PList_remItem(PList * lst, MDB_MeshEntity * ent)
  {
    std::vector<MDB_MeshEntity *>::iterator eIter = lst->entities.begin();
    for (; eIter != lst->entities.end() ; eIter++) {
      if ( *eIter == ent ) lst->entities.erase(eIter);
    }
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

}
