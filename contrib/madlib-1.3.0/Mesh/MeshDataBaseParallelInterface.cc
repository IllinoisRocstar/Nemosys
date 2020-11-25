// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: J.-F. Remacle, C. Dobrzynski, K. Hillewaert, G. Compere
// -------------------------------------------------------------------

#include "assert.h"
#include "MeshDataBase.h"
#include "MeshDataBaseIO.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseParallelInterface.h"
#include "MeshDataBaseParallelIO.h"
#include "MeshDataBaseCommCheck.h"
#include "MeshDataBaseComm.h"
#include "MeshDataBaseCommPeriodic.h"
#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#endif

#include <set>
#include <iterator>
#include <algorithm>

namespace MAd {

  struct point_comm
  {
    pVertex p;     // local pointer
    int nID;       // myrank
    int numGlobal; // Id global for the vertex
  };
  
  struct edge_comm {
    int distProc;  // remote processor 
    pEdge pe;      // edge 
    pVertex p1,p2; // vertices on remote processor
    int id1,id2;   // local identity tags
  };  
  
  struct face_comm
  {
    int     distProc;        // remote processor 
    pFace   pf;              // face
    pVertex p1,p2,p3,p4;     // vertices on remote processor
    int     id1,id2,id3,id4; // local identity tags 
  };

  // -------------------------------------------------------------------
  /*! \brief write mesh in parallel format \ingroup internal */
#ifdef PARALLEL
  void M_writeParallel(pMesh m, const char * filename, int version)
  {
    SaveGmshMeshParallel(m, filename, version);
  }
#endif

  // -------------------------------------------------------------------
  /*! \brief flags true if entity pe is on a partition boundary \ingroup internal
    \warning only valid after creation of interface exchange information (based
    on data with tag "RemotePoint") */
  
  bool EN_isInterface(pEntity pe)
  {
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    void *temp_ptr;
    int isInterface = EN_getDataPtr(pe, tagData, &temp_ptr);
    if (isInterface) return true;
    return false;
  }
  
  // -------------------------------------------------------------------
  /*! \brief flags true if vertex pv is on a partition boundary \ingroup parallel
    \warning only valid after V_createInfoInterface (based on EN_isInterface) */
  
  bool V_isInterface(pVertex pv) { return EN_isInterface((pEntity) pv);}
  
  // -------------------------------------------------------------------
  /*! \brief flags true if edge pe is on a partition boundary \ingroup parallel
    \warning only valid after E_createInfoInterface (based on EN_isInterface)  */
  
  bool E_isInterface(pEdge pe)   { return EN_isInterface((pEntity) pe);}
  
  // -------------------------------------------------------------------
  /*! \brief flags true if face pf is on a partition boundary \ingroup parallel
    \warning only valid after F_createInfoInterface  (based on EN_isInterface) */

  bool F_isInterface(pFace pf)   { return EN_isInterface((pEntity) pf);}

  // -------------------------------------------------------------------
  /*!  \brief Lists all processors on which the vertex pv has a copy \ingroup internal 
    The data with tag "RemotePoint" is created during
    V_createInfoInterface and stores a list of pairs consisting of processor
    id and vertex pointer
  */
  
  int V_listInterface(pVertex pv, std::vector<int>* distProcs)
  {
    pMeshDataId tagVertex = MD_lookupMeshDataId("RemotePoint");
    
    void *temp_ptr;
    int isInterface = EN_getDataPtr((pEntity) pv, tagVertex, &temp_ptr);
    if(isInterface) {
      const std::vector<std::pair<int,pVertex> > *recup = 
        (const std::vector<std::pair<int,pVertex> > *) temp_ptr;
      
      std::set<int> parts;
      std::vector<std::pair<int,pVertex> >::const_iterator pIter = recup->begin();
      for (;pIter!=recup->end();++pIter) parts.insert(pIter->first);
      
      distProcs->insert(distProcs->begin(),parts.begin(),parts.end());
      
    }
    return distProcs->size();
  }
  
  // -------------------------------------------------------------------
  /*!  \brief Returns true if the node is a (periodic) copy of the node with global id rId \ingroup internal */
  /*! This information is based on data read with the mesh and hence will always be correct */
  
  bool V_corresponds(pVertex pv,int rId) 
  {
    if (EN_id((pEntity) pv) == rId) return true;
    
    pMeshDataId tagVertex = MD_lookupMeshDataId("PeriodicVertexID");
    void* tmpptr = NULL;
    if (!EN_getDataPtr((pEntity) pv,tagVertex,&tmpptr)) return false;
    std::map<int,std::vector<int> >* conn = (std::map<int,std::vector<int> >*) tmpptr;
    return (conn->find(rId) != conn->end());
  }

  // -------------------------------------------------------------------
  
  /*! \brief return true if the edge corresponds (possibly) \ingroup internal */
  /*! \warning error for edges transformed onto oneself */
  
  bool E_corresponds(pEdge pe,int rv1,int rv2) {

    pVertex pv1 = E_vertex(pe,0);
    pVertex pv2 = E_vertex(pe,1);


    int lv1 = EN_id((pEntity) pv1);
    int lv2 = EN_id((pEntity) pv2);
    
    int lmin = std::min(lv1,lv2);
    int lmax = std::max(lv1,lv2);

    int rmin = std::min(rv1,rv2);
    int rmax = std::max(rv1,rv2);

    // identical edge - most of the cases
    
    if (lmin == rmin && lmax == rmax) return true;

    // only one shared node -> problems
    
    if (lmin == rmax && lmax != rmax) return false;
    if (lmin != rmin && lmax == rmax) return false;
    
    // all other cases should be ok ? (Koen)
    
    if (V_corresponds(pv1,rv1) && V_corresponds(pv2,rv2)) return true;
    if (V_corresponds(pv1,rv2) && V_corresponds(pv2,rv1)) return true;
    
  }

  // -------------------------------------------------------------------
  /*! \brief return true if the face corresponds \ingroup internal */
  /*! \warning not yet fully functional outside of F_createInfoInterface */
  /*! <ul>
      <li> we do not really verify periodic nodes 
      <li> error for edges transformed onto oneself
      </ul>
  */  

  bool F_corresponds(pFace pf,int rv1,int rv2,int rv3,int rv4) {
    
    std::set<int> lverts;

    for (int i=0;i<F_numVertices(pf);i++) {
      lverts.insert(EN_id((pEntity) F_vertex(pf,i)));
    }
    
    std::set<int> rverts;
    rverts.insert(rv1);
    rverts.insert(rv2);
    rverts.insert(rv3);
    if (rv4 != -1) rverts.insert(rv4);
    
    // obvious checks ... 
    if (rverts.size() != lverts.size()) return false;

    // same face ... 
    if (rverts == lverts) return true;

    // only some shared nodes is not ok ? 
    std::set<int> nonShared;
    std::set_difference(lverts.begin(),lverts.end(),
                        rverts.begin(),rverts.end(),
                        std::insert_iterator< std::set<int> > (nonShared,
                                                               nonShared.begin()));
    
    if (nonShared.size() != lverts.size()) return false;
    return true;
  }

  // -------------------------------------------------------------------
  /*!  \brief Returns true if the node is a (periodic) copy of the node with global id rId \ingroup parallel 
    
    This information is based on data read with the mesh and hence will always
    be correct; only valid for check on local partition */
  
  bool V_corresponds(pVertex p1,pVertex p2) 
  {
    return (V_corresponds(p1,EN_id((pEntity) p2)) && 
            V_corresponds(p2,EN_id((pEntity) p1)));
  }

  // --------------------------------------------------------------------

  /*! \brief returns true if the node is flagged as periodic \ingroup parallel

    This information is based on data read with the mesh and hence will always
    be correct */

  bool V_isPeriodic(pVertex pv) {
    
    pMeshDataId tagVertex = MD_lookupMeshDataId("PeriodicVertexID");
    void* tmpptr = NULL;
    if (!EN_getDataPtr((pEntity) pv,tagVertex,&tmpptr)) return false;
    return true;
  }
    
  // -------------------------------------------------------------------
  /*!\brief  An edge flagged as potentially on a parallel interface  \ingroup internal

    \warning the check is not correct for 2D meshes with internal edges */    
  
  bool E_isPotentialInterface(pMesh mesh, pEdge pe)
  {

    switch (M_dim(mesh)) {
    case 2:
      {     
      
        // --- check that the edge has less than two neighbouring faces ---
        
        int numf = pe->numfaces();
        if( numf == 2 ) return false;
        
        // --- check that the edge is not alone or
        //     with multiple neighbouring faces (debug) ---
        
        assert(numf == 1);
        return true;
        
        // KH: for the moment an edge is always classified on an edge ... (nullmodel)
        //     hence the following check is not usable ... 
        /// 
        // if edge is classified on an edge, still two possibilities
        // 
        //   1. periodic edge: check whether or not nodes are periodic - ok
        //   2. internal edge: not yet treated correctly
        
        
        if ( E_whatInType(pe) == 1 ) {
          for (int i=0;i<2;i++) if (!V_isPeriodic(E_vertex(pe,i))) return false;
          return true;
        }
        
      
        break;
      }
    case 3:
      {
        // --- check that the edge has at least one neighbouring face that is potentially shared ---
        
        int numf = pe->numfaces();
        for (int iF=0; iF < numf; iF++) {
          if ( F_isPotentialInterface(mesh,E_face(pe,iF))) return true;
        }
        break;
      }
    }

    return false;
  }
  
  // -------------------------------------------------------------------

//   bool E_isInterface(pMesh mesh, pEdge pe, int * distProc, 
//                      std::vector<pVertex>* distVt)
//   {
//     if (distProc) *distProc = -1;
//     if (distVt)   (*distVt).clear();

//     if ( !E_isPotentialInterface(mesh,pe) ) return false;

//     // GCREMARK: the rest of the function gives informations about the 
//     //           distant partition(s) but is just a cross-check.

//     // --- check that the two nodes are on a parallel boundary ---
//     // KHREMARK :
//     // - this is not conclusive, we should verify that the edge exists on the other partitions
//     // - on the other hand, if the edge is possibly on more than two partitions,
//     //   we will not flag all remote partitions
  
//     pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
//     pVertex pv1 = pe->p1;
//     pVertex pv2 = pe->p2;
//     void * temp_ptr1, *temp_ptr2;
//     int isInterface1 = EN_getDataPtr((pEntity) pv1 , tagData, &temp_ptr1);
//     int isInterface2 = EN_getDataPtr((pEntity) pv2 , tagData, &temp_ptr2);
//     if( !(isInterface1 && isInterface2) ) {
//       Msg(MDB_FATAL,"Error: found an edge which nodes are not both on a parallel boundary although it fits the conditions to be a parallel boundary edge\n");
//     }

//     // --- check that the two nodes are on the same parallel boundary ---
//     const std::vector<std::pair<int,pVertex> > *recup1 = 
//       (const std::vector<std::pair<int,pVertex> > *) temp_ptr1;
//     const std::vector<std::pair<int,pVertex> > *recup2 = 
//       (const std::vector<std::pair<int,pVertex> > *) temp_ptr2;
        
//     std::vector<std::pair<int,pVertex> >::const_iterator rIter1 = recup1->begin();
//     std::vector<std::pair<int,pVertex> >::const_iterator rIter2 = recup2->begin();
        
//     for (; rIter1 != recup1->end(); rIter1++) {
//       for (; rIter2 != recup2->end(); rIter2++) {
//         if ( (*rIter1).first == (*rIter2).first ) {
//           if (distProc) *distProc = (*rIter1).first;
//           if (distVt) {
//             (*distVt).clear();
//             (*distVt).push_back((*rIter1).second);
//             (*distVt).push_back((*rIter2).second);
//           }
//           return true;
//         }
//       }
//       rIter2 = recup2->begin();
//     }

//     Msg(MDB_FATAL,"Error: found an edge which nodes have no common parallel boundary although it fits the conditions to be a parallel boundary edge\n");

//     return false;
//   }


  // ---------------------------------------------------------------------------
  /*!\brief Verifies whether edge pe is located at a parallel interface boundary 
    and provides information about all possible copies \ingroup internal */
  
  bool E_isPotentialInterface(pMesh mesh, pEdge pe, std::vector<edge_comm>& distVt)

  {
    distVt.clear();
    
    if ( !E_isPotentialInterface(mesh,pe) ) return false;

    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    pVertex pv1 = pe->p1;
    pVertex pv2 = pe->p2;
    void * temp_ptr1, *temp_ptr2;
    int isInterface1 = EN_getDataPtr((pEntity) pv1 , tagData, &temp_ptr1);
    int isInterface2 = EN_getDataPtr((pEntity) pv2 , tagData, &temp_ptr2);
  
    if( !(isInterface1 && isInterface2) ) return false;

    // --- check that the two nodes are on the same parallel boundary ---
    const std::vector<std::pair<int,pVertex> > *recup1 = 
      (const std::vector<std::pair<int,pVertex> > *) temp_ptr1;
    const std::vector<std::pair<int,pVertex> > *recup2 = 
      (const std::vector<std::pair<int,pVertex> > *) temp_ptr2;

    std::vector<std::pair<int,pVertex> >::const_iterator rIter1 = recup1->begin();
    std::vector<std::pair<int,pVertex> >::const_iterator rIter2 = recup2->begin();
  
    for (; rIter1 != recup1->end(); rIter1++) {
      for (; rIter2 != recup2->end(); rIter2++) {  
        if ( (*rIter1).first == (*rIter2).first ) {
          edge_comm ee;
          ee.distProc = (*rIter1).first;
          ee.p1 = (*rIter1).second;
          ee.p2 = (*rIter2).second;
          ee.id1 = EN_id((pEntity) E_vertex(pe,0));
          ee.id2 = EN_id((pEntity) E_vertex(pe,1));
          ee.pe = pe;
          distVt.push_back(ee);
        }
      }
      rIter2 = recup2->begin();
    }
    return (!distVt.empty());
  }

  // -------------------------------------------------------------------
  /*! \brief Verify whether pf may potentially be an interface \ingroup internal 
    \warning will flag false for internal boundaries !!! */ 
  
  // -------------------------------------------------------------------

  bool F_isPotentialInterface(pMesh mesh, pFace pf)
  {
    // bool check = true;
    //     for (size_t i=0;i<F_numVertices(pf);i++) check = check && V_isInterface(F_vertex(pf,i));
    //     return check;

    // --- check that the mesh is 3D ---
    if ( M_dim(mesh) != 3 ) return false;
  
    // --- check that the face has less than two neighbouring regions ---

    int numtet = pf->getNbRegions();
    if (numtet == 2) return false;
    return true;

    // OK : for the moment a face is always classified on a geometric face (nullmodel ...)
    
    // --- if the face is classified on a topological face, we can still have 
    //     1. a periodic face: all nodes are periodic - ok
    //     2. an internal boundary - check is not ok 
    
    int modelDim = GEN_type(pf->g);

    if (modelDim == 2) {
      for (int i=0;i<F_numVertices(pf);i++) {
        if (!V_isPeriodic(F_vertex(pf,i))) return false;
      }
      return true;
    }
    
    // --- check that the face is not alone (debug) ---
    assert(numtet);
    return true;
  }


  // -------------------------------------------------------------------
#ifdef PARALLEL
  // bool F_isInterface(pMesh mesh, pFace pface, int * distProc,
  //                    std::vector<pVertex>* distVt)
  // {
  //   if (distProc) *distProc = -1;
  //   if (distVt)   (*distVt).clear();

  //   if ( !F_isInterface(mesh,pface) ) return false;

  //   // GCREMARK: the rest of the function gives informations about the 
  //   //           distant partition but is just a cross-check.

  //   // --- check that the three nodes are on a parallel boundary ---
  //   pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
  //   pVertex nod[4];


  //   switch (pface->getNbNodes()) {

  //   case 3:
  //     {
  //       pface->getNodes(nod);
  //       void * temp_ptr1, *temp_ptr2, *temp_ptr3;
  //       int isInterface1 = EN_getDataPtr((pEntity) nod[0], tagData, &temp_ptr1);
  //       int isInterface2 = EN_getDataPtr((pEntity) nod[1], tagData, &temp_ptr2);
  //       int isInterface3 = EN_getDataPtr((pEntity) nod[2], tagData, &temp_ptr3);
      
      
  //       if( !(isInterface1 && isInterface2 && isInterface3) ) {
  //         Msg(MDB_FATAL,"Error: found a face which nodes are not all on a parallel boundary although it fits the conditions to be a parallel boundary face\n");
  //       }
      
  //       // --- check that the three nodes are on the same parallel boundary ---
  //       const std::vector<std::pair<int,pVertex> > *recup1 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr1;
  //       const std::vector<std::pair<int,pVertex> > *recup2 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr2;
  //       const std::vector<std::pair<int,pVertex> > *recup3 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr3;
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter1 = recup1->begin();
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter2 = recup2->begin();
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter3 = recup3->begin();
  //       for (; rIter1 != recup1->end(); rIter1++) {
  //         for (; rIter2 != recup2->end(); rIter2++) {
  //           for (; rIter3 != recup3->end(); rIter3++) {
  //             if ( (*rIter1).first == (*rIter2).first &&
  //                  (*rIter2).first == (*rIter3).first ) {
  //               if (distProc) *distProc = (*rIter1).first;
  //               if (distVt) {
  //                 (*distVt).push_back((*rIter1).second);
  //                 (*distVt).push_back((*rIter2).second);
  //                 (*distVt).push_back((*rIter3).second);
  //               }
  //               return true;
  //             }
  //           }
  //           rIter3 = recup3->begin();
  //         }
  //         rIter2 = recup2->begin();
  //       }
  //       break;
  //     }
    
  //   case 4:
  //     {
  //       pface->getNodes(nod);
  //       void * temp_ptr1, *temp_ptr2, *temp_ptr3, *temp_ptr4;
  //       int isInterface1 = EN_getDataPtr((pEntity) nod[0], tagData, &temp_ptr1);
  //       int isInterface2 = EN_getDataPtr((pEntity) nod[1], tagData, &temp_ptr2);
  //       int isInterface3 = EN_getDataPtr((pEntity) nod[2], tagData, &temp_ptr3);
  //       int isInterface4 = EN_getDataPtr((pEntity) nod[3], tagData, &temp_ptr4);
      
    
  //       if( !(isInterface1 && isInterface2 && isInterface3 && isInterface4) ) {
  //         Msg(MDB_FATAL,"Error: found a face which nodes are not all on a parallel boundary although it fits the conditions to be a parallel boundary face\n");
  //       }
    
  //       // --- check that the three nodes are on the same parallel boundary ---
  //       const std::vector<std::pair<int,pVertex> > *recup1 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr1;
  //       const std::vector<std::pair<int,pVertex> > *recup2 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr2;
  //       const std::vector<std::pair<int,pVertex> > *recup3 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr3;
  //       const std::vector<std::pair<int,pVertex> > *recup4 = 
  //         (const std::vector<std::pair<int,pVertex> > *) temp_ptr4;
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter1 = recup1->begin();
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter2 = recup2->begin();
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter3 = recup3->begin();
  //       std::vector<std::pair<int,pVertex> >::const_iterator rIter4 = recup4->begin();
  //       for (; rIter1 != recup1->end(); rIter1++) {
  //         for (; rIter2 != recup2->end(); rIter2++) {
  //           for (; rIter3 != recup3->end(); rIter3++) {
  //             for (; rIter4 != recup4->end(); rIter4++) {
  //               if ( (*rIter1).first == (*rIter2).first &&
  //                    (*rIter2).first == (*rIter3).first &&
  //                    (*rIter3).first == (*rIter4).first ) {
  //                 if (distProc) *distProc = (*rIter4).first;
  //                 if (distVt) {
  //                   (*distVt).push_back((*rIter1).second);
  //                   (*distVt).push_back((*rIter2).second);
  //                   (*distVt).push_back((*rIter3).second);
  //                   (*distVt).push_back((*rIter4).second);
  //                 }
  //                 return true;
  //               }
  //             }
  //             rIter4 = recup4->begin();
  //           }
  //           rIter3 = recup3->begin();
  //         }
  //         rIter2 = recup2->begin();
  //       }
  //       break;
  //     }
  //   }
    
  //   Msg(MDB_FATAL,"Error: found a face which nodes have no common parallel boundary although it fits the conditions to be a parallel boundary face\n");

  //   return false;
  // }

#endif

  /*! \brief flag face pface as potential interface and list possible communications \ingroup internal */

  bool F_isPotentialInterface(pMesh mesh, pFace pface, std::vector<face_comm>& distFace)
  {

    if ( !F_isPotentialInterface(mesh,pface) ) return false;

    // GCREMARK: the rest of the function gives informations about the 
    //           distant partition but is just a cross-check.

    // --- check that the three nodes are on a parallel boundary ---
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    pVertex nod[4];

    distFace.clear();


    switch (pface->getNbNodes()) {

    case 3:
      {
        pface->getNodes(nod);
        void * temp_ptr1, *temp_ptr2, *temp_ptr3;
        int isInterface1 = EN_getDataPtr((pEntity) nod[0], tagData, &temp_ptr1);
        int isInterface2 = EN_getDataPtr((pEntity) nod[1], tagData, &temp_ptr2);
        int isInterface3 = EN_getDataPtr((pEntity) nod[2], tagData, &temp_ptr3);
      
        int id1 = EN_id((pEntity) nod[0]);
        int id2 = EN_id((pEntity) nod[1]);
        int id3 = EN_id((pEntity) nod[2]);

        if( !(isInterface1 && isInterface2 && isInterface3) ) {
          return false;
          //Msg(MDB_FATAL,"Error: found a face which nodes are not all on a parallel boundary although it fits the conditions to be a parallel boundary face\n");
        }
      
        // --- check that the three nodes are on the same parallel boundary ---
        const std::vector<std::pair<int,pVertex> > *recup1 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr1;
        const std::vector<std::pair<int,pVertex> > *recup2 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr2;
        const std::vector<std::pair<int,pVertex> > *recup3 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr3;
        std::vector<std::pair<int,pVertex> >::const_iterator rIter1 = recup1->begin();
        std::vector<std::pair<int,pVertex> >::const_iterator rIter2 = recup2->begin();
        std::vector<std::pair<int,pVertex> >::const_iterator rIter3 = recup3->begin();
        for (; rIter1 != recup1->end(); rIter1++) {
          for (; rIter2 != recup2->end(); rIter2++) {
            for (; rIter3 != recup3->end(); rIter3++) {
              if ( (*rIter1).first == (*rIter2).first &&
                   (*rIter2).first == (*rIter3).first ) {

                face_comm fc;
              
                fc.distProc = (*rIter1).first;

                fc.p1 = (*rIter1).second;
                fc.p2 = (*rIter2).second;
                fc.p3 = (*rIter3).second;
                fc.p4 = NULL;
                
                fc.id1 = id1;
                fc.id2 = id2;
                fc.id3 = id3;
                fc.id4 = -1;

                fc.pf = pface;
              
                distFace.push_back(fc);
              }
            }
            rIter3 = recup3->begin();
          }
          rIter2 = recup2->begin();
        }
        break;
      }
    
    case 4:
      {
        pface->getNodes(nod);
        void * temp_ptr1, *temp_ptr2, *temp_ptr3, *temp_ptr4;
        int isInterface1 = EN_getDataPtr((pEntity) nod[0], tagData, &temp_ptr1);
        int isInterface2 = EN_getDataPtr((pEntity) nod[1], tagData, &temp_ptr2);
        int isInterface3 = EN_getDataPtr((pEntity) nod[2], tagData, &temp_ptr3);
        int isInterface4 = EN_getDataPtr((pEntity) nod[3], tagData, &temp_ptr4);
        
        int id1 = EN_id((pEntity) nod[0]);
        int id2 = EN_id((pEntity) nod[1]);
        int id3 = EN_id((pEntity) nod[2]);
        int id4 = EN_id((pEntity) nod[3]);
    
        if( !(isInterface1 && isInterface2 && isInterface3 && isInterface4) ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "found a face which nodes are not all on a parallel boundary although it fits the conditions to be a parallel boundary face");
        }
    
        // --- check that the three nodes are on the same parallel boundary ---
        const std::vector<std::pair<int,pVertex> > *recup1 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr1;
        const std::vector<std::pair<int,pVertex> > *recup2 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr2;
        const std::vector<std::pair<int,pVertex> > *recup3 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr3;
        const std::vector<std::pair<int,pVertex> > *recup4 = 
          (const std::vector<std::pair<int,pVertex> > *) temp_ptr4;
        std::vector<std::pair<int,pVertex> >::const_iterator rIter1 = recup1->begin();
        std::vector<std::pair<int,pVertex> >::const_iterator rIter2 = recup2->begin();
        std::vector<std::pair<int,pVertex> >::const_iterator rIter3 = recup3->begin();
        std::vector<std::pair<int,pVertex> >::const_iterator rIter4 = recup4->begin();
        for (; rIter1 != recup1->end(); rIter1++) {
          for (; rIter2 != recup2->end(); rIter2++) {
            for (; rIter3 != recup3->end(); rIter3++) {
              for (; rIter4 != recup4->end(); rIter4++) {
                if ( (*rIter1).first == (*rIter2).first &&
                     (*rIter2).first == (*rIter3).first &&
                     (*rIter3).first == (*rIter4).first ) {
                
                  face_comm fc;
                
                  fc.distProc = (*rIter1).first;

                  fc.p1 = (*rIter1).second;
                  fc.p2 = (*rIter1).second;
                  fc.p3 = (*rIter2).second;
                  fc.p4 = (*rIter3).second;
                  
                  fc.id1 = id1;
                  fc.id2 = id2;
                  fc.id3 = id3;
                  fc.id4 = id4;

                  fc.pf = pface;
                
                  distFace.push_back(fc);
                
                }
              }
              rIter4 = recup4->begin();
            }
            rIter3 = recup3->begin();
          }
          rIter2 = recup2->begin();
        }
        break;
      }
    }
    
    // Msg(MDB_FATAL,"Error: found a face which nodes have no common parallel boundary although it fits the conditions to be a parallel boundary face\n");

    return (distFace.size() != 0);
  }

  // -------------------------------------------------------------------
  /*! \brief Flag region if located at an interface. \ingroup parallel
    \warning As we partition along edges in 2D and faces in 3D, the answer is
    always no  */

  bool R_isInterface(pRegion pr) {
    return false;
  }

  // -------------------------------------------------------------------
  // do a tentative connection - since we send all potential interface nodes
  // to all processors, we send the minimum information : the node id 
  // 
  // result is that we store the connected partitions for a shared node
  // more detailed communication follows 
  
  
  void createRemotePointLists (pMesh mesh, pMeshDataId tagVertex)
  {
#ifdef PARALLEL
    int myrank = 0;
    int nbproc = 1;
#endif

    std::set<pVertex> bdryNodes;
    std::set<int>     bdryNodeId;

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // --- list local nodes on parallel boundaries ---
   int meshDim = M_dim(mesh);

    if ( meshDim == 3 ) {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        if ( F_isPotentialInterface(mesh,pf) ) {
          for (int iv = 0; iv < F_numVertices(pf); iv++) {
            pVertex pv = F_vertex(pf,iv);
            bdryNodes .insert(pv);
            bdryNodeId.insert(pv->iD);
          }
        }
      }
    } else {
      EIter eit = M_edgeIter(mesh);
      while ( pEdge pe = EIter_next(eit) ) {
        if ( E_isPotentialInterface(mesh,pe) ) {
          for (int iv = 0; iv < 2; iv++) {
            pVertex pv = E_vertex(pe,iv);            
            bdryNodes.insert(pv);
            bdryNodeId.insert(pv->iD);
          }
        }
      }
    }
#endif

    /* KH : temporary ? addition of all periodic nodes */ 
    
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicVertexID");  
    
    VIter vit = M_vertexIter(mesh);
    
    while (pVertex pv = VIter_next(vit)) {
      
      void* tmp_periodic = NULL;
      int isPeriodic = EN_getDataPtr((pEntity) pv,tagPeriodic,&tmp_periodic);

      if (isPeriodic) {
      
        void* tmp_remote = NULL;
        int hasRemote  = EN_getDataPtr((pEntity) pv,tagVertex,&tmp_remote);
        
        if (!hasRemote) {
          std::vector<std::pair<int,pVertex> >* remote = new std::vector<std::pair<int,pVertex> >;
          EN_attachDataPtr((pEntity) pv,tagVertex,remote);
        }
        else {
          std::vector<std::pair<int,pVertex> >* remote = (std::vector<std::pair<int,pVertex> >*) tmp_remote;
          remote->clear();
        }
        
        std::map<int,std::vector<int> > * conn = (std::map<int,std::vector<int> > *) tmp_periodic;
        std::map<int,std::vector<int> >::iterator citer = conn->begin();
        bdryNodes.insert(pv);
        for (;citer != conn->end();++citer) bdryNodeId.insert(citer->first);
      }
    }
    
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    
    // --- send list ---
    
    int* sendcounts = new int[nbproc];
    int size = bdryNodeId.size()*sizeof(int);
    
    void* msg = malloc(size);
    int * bufCast = reinterpret_cast<int *>(msg);
    std::set<int>::const_iterator cIter = bdryNodeId.begin();
    std::set<int>::const_iterator cLast = bdryNodeId.end();
    for (;cIter!=cLast;++cIter) *(bufCast++) = *cIter;
    
    for (int iProcDest=0; iProcDest < nbproc; iProcDest++) {
      if ( iProcDest == myrank ) {
        sendcounts[iProcDest] = 0;
        continue;    
      }
      else {
        sendcounts[iProcDest] = 1;
        void *buf = AP_alloc(iProcDest,669,size);
        memcpy(buf,msg,size);
        AP_send(buf);
      }
    }

    free(msg);

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
    AP_flush();

    // --- receive lists ---

    int message = 0;
    int count;
    while (!AP_recv_count(&count) || message < count) {

      void *msg;
      int from;
      int tag;
      int size;
      int rc;
    
      rc = AP_recv (MPI_ANY_SOURCE, 669, AP_BLOCKING|AP_DROPOUT,
                    &msg, &size, &from, &tag);
    
      if (rc) {
    
        int* msgCast = reinterpret_cast<int*>(msg);

        std::set<int> distNodes;
        int nbDistNodes = size / sizeof(int);
        for (int iDistN = 0; iDistN < nbDistNodes; iDistN++) {
          distNodes.insert(msgCast[iDistN]);
        }
        
        std::set<pVertex>::const_iterator vIter = bdryNodes.begin();
        std::set<pVertex>::const_iterator vLast = bdryNodes.end();
        
        for (;vIter!=vLast;++vIter) {

          pVertex pv = *vIter;
          int     id = EN_id((pEntity) pv);

          if (distNodes.find(id) != distNodes.end()) {
            void *temp_ptr;
            int exist = EN_getDataPtr((pEntity) pv, tagVertex, &temp_ptr);
            std::vector<std::pair<int,pVertex> > * parts = NULL;
            if (exist) {
              parts = reinterpret_cast<std::vector<std::pair<int,pVertex> >*>(temp_ptr);
            } else if (!exist) {
              parts = new std::vector<std::pair<int,pVertex> >;
              EN_attachDataPtr((pEntity) pv, tagVertex, parts);
            }
            parts->push_back(std::pair <int, MDB_Point*>(from,NULL));
          }
        } 
        message++;
        AP_free(msg);
      }
    }
    
    // --- wait until all have finished sending ---
  
    AP_check_sends(AP_WAITALL);
    delete [] sendcounts;
    
    // GCDEBUG
    //   MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  // -------------------------------------------------------------------
  
  void V_createInfoInterface(pMesh mesh, pMeshDataId tagVertex)
  {
        
    int mysize = 1;
    int myrank = 0;

    size_t nbConnections = 0;
    
#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

    // --- Delete previous information ---
    VIter vit = M_vertexIter(mesh);
    while ( pVertex pv = VIter_next(vit) ) {
      void * tmp_ptr;
      if( EN_getDataPtr((pEntity) pv, tagVertex, &tmp_ptr) ) {
        EN_deleteData((pEntity) pv, tagVertex);
      }
    }
    VIter_delete(vit);

    // --- Create lists of distant procs ---

    createRemotePointLists(mesh, tagVertex);

    // --- Send node pointers ---
    int *sendcounts = new int[mysize];
    for(int i=0; i<mysize; i++) sendcounts[i] = 0;
    
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicVertexID");
    //pMeshDataId tagTransfo  = MD_lookupMeshDataId("PeriodicTransformation");
  
    vit = M_vertexIter(mesh);
    while ( pVertex pv = VIter_next(vit) ) {

      std::vector<int> distProcTab;
      //int nbDistProc = V_listInterface(pv, &distProcTab);

      std::set<int> distProcs;
      distProcs.insert(distProcTab.begin(),distProcTab.end());
      
      void* tmp_remote = NULL;
      int haveRemote = EN_getDataPtr((pEntity) pv,tagVertex,&tmp_remote);
      std::vector<std::pair<int,pVertex> >* remote = NULL;
      if (haveRemote) {
        remote  = (std::vector<std::pair<int,pVertex> >*) (tmp_remote);
        remote->clear();
      }
      
      void* tmp_periodic = NULL;
      int havePeriodic = EN_getDataPtr((pEntity) pv,tagPeriodic,&tmp_periodic);
      std::map<int,std::vector<int> >* periodic = NULL;

      if (havePeriodic) {
        
        if (!haveRemote) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Periodic node %d without connection table",
                                      EN_id((pEntity) pv));
        }
        
        // allocate and attach transformation table 
        
        //         int haveTransfo = EN_getDataPtr((pEntity) cvtx, tagTransfo, &tmp_transfo);
        //         std::map<std::pair<int,pEntity>,std::vector<int> >* trafo = NULL; 
        //         if (haveTransfo) {
        //           trafo = (std::map<std::pair<int,pEntity>,std::vector<int> > *) tmp_transfo;
        //           trafo->clear();
        //         } 
        //         else {
        //           trafo = new std::map<std::pair<int,pEntity>,std::vector<int> >;
        //           EN_attachDataPtr((pEntity) cvtx,periodicTransfoTag,trafo);
        //         }
        
        // insert data on local processor
        
        periodic = (std::map<int,std::vector<int> >*) tmp_periodic;
        std::map<int,std::vector<int> >::const_iterator citer = periodic->begin();
        for (;citer!=periodic->end();++citer) {
          pVertex cvtx = mesh->find_point(citer->first);
          if (cvtx && cvtx != pv) {
            remote->push_back(std::make_pair(myrank,cvtx));
            nbConnections++;
            //     trafo->insert(std::make_pair(std::make_pair(myRank,cvtx),citer->second));
          }
        }
      }
      
#ifdef PARALLEL
      
      int nbConnections = 1;
      if (periodic) nbConnections += periodic->size();
      int msgSize = sizeof(pVertex) + (3 + nbConnections ) * sizeof(int);
      void* tmp = malloc(msgSize);
      
      pVertex *vbuf = (pVertex *) tmp;

      int lId = EN_id((pEntity) pv);
      
      *(vbuf++) = pv;
      int* iv = (int*) (vbuf);
      
      *(iv++) = GEN_tag (V_whatIn(pv));
      *(iv++) = GEN_type(V_whatIn(pv));
      *(iv++) = nbConnections;
      *(iv++) = lId;
      
      if (periodic) {
        std::map<int,std::vector<int> >::const_iterator cIter = periodic->begin();
        for (;cIter!=periodic->end();++cIter) {
          int rId = cIter->first;
          *(iv++) = (rId == lId) ? -rId : rId; // avoid sending once again to oneself
        }
      }
      
      // make sure that we send only once to different processors

      std::set<int>::iterator dIter = distProcs.begin();
      for (;dIter!=distProcs.end();++dIter) {

        int distProc = *dIter;
        if (distProc != myrank) {    
          void *buf = AP_alloc(distProc,445,msgSize);
          memcpy(buf,tmp,msgSize);
          AP_send(buf);
          sendcounts[distProc]++;
        }
      }
      
      free(tmp);
#endif
    }
    VIter_delete(vit);
    
#ifdef PARALLEL
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
    AP_flush();

    // --- Receive node pointers ---
    int message=0, count;
    while (!AP_recv_count(&count) || message<count) {
      
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      
      recv = AP_recv(MPI_ANY_SOURCE, 445, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        
        message++;
        pVertex* vbuf = reinterpret_cast<pVertex*>(msg);
        pVertex precv = *(vbuf++);
        int* ibuf = reinterpret_cast<int*>(vbuf);
        
        int gTag   = *(ibuf++);
        int gDim   = *(ibuf++);
        int nbConn = *(ibuf++);
        
        
        for (int iConn=0;iConn<nbConn;iConn++) {
          
          int pId = *(ibuf++);
          pVertex vtx = mesh->find_point(pId);
          
          if (vtx) {     
            
            void* tmp_remote =  NULL;
            int isInterface = EN_getDataPtr((pEntity) vtx, tagVertex, &tmp_remote);
            assert(isInterface);
            std::vector<std::pair<int,pVertex> > *remote = (std::vector<std::pair<int,pVertex> >*) tmp_remote;  
            remote->push_back(std::pair<int,pVertex>(from,precv));
            
          }
        }
        AP_free(msg);
      }
    }
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;

#endif
    
    MDB_CommCheck cc;
    exchangeDataOnVertices(mesh,cc);

  }
  
//   // -------------------------------------------------------------------
// #ifdef PARALLEL
//   void MDB_Mesh::bdryLinkRemotePoint()
//   {
//     int mysize,myrank;
//     MPI_Comm_size(MPI_COMM_WORLD, &mysize);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);   

//     /*send pointeurs*/
//     int *sendcounts = new int[mysize];
//     for(int i=0;i<mysize;i++)sendcounts[i]=0;
  
//     pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
//     VIter vit = M_vertexIter(this);
//     pVertex pv;  
//     int nsend=0;
//     while ((pv = VIter_next(vit))) {
//       void *temp_ptr; 
//       int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
    
//       if(isInterface) {
//         std::vector<std::pair<int , MDB_Point*> > *recup = (std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
//         int numGlobal = EN_id((pEntity)pv);
//         for(unsigned int j=0 ; j<(*recup).size() ; j++) {
//           nsend++;
//           int remoteID  = (*recup)[j].first;
//           assert(remoteID != myrank);
//           void *buf = AP_alloc(remoteID,444,sizeof(point_comm));
//           point_comm *castbuf = (point_comm *) buf;
//           castbuf->p         = pv;
//           castbuf->nID       = myrank;
//           castbuf->numGlobal = numGlobal;
//           AP_send(buf);
//           sendcounts[remoteID]++;   
//         }
//       }
//     }	
//     VIter_delete(vit);   

//     AP_flush();
//     AP_check_sends(AP_NOFLAGS);
//     // AP_reduce_nsends(sendcounts);

//     /*receive pointers*/
//     int message=0;
//     // while (!AP_recv_count(&count) || message<count) {
//     while (message<nsend) {
//       void *msg;
//       int  from;
//       int  tag;
//       int  size;
//       int  recv;
//       recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
//                      &msg, &size, &from, &tag);
//       if(recv) {
//         message++;
//         point_comm * castbuf = (point_comm*) msg;
//         pVertex precv (0);
//         precv = castbuf -> p;
//         assert(precv);
//         int recvID = castbuf -> nID;
//         int numrecv = castbuf -> numGlobal;
      
//         pVertex p = this->find_point(numrecv);
//         if(p) {
//           void *temp_ptr; 
//           int isInterface = EN_getDataPtr((pEntity) p , tagData, &temp_ptr);
//           assert(isInterface);
//           std::vector<std::pair<int , MDB_Point*> > *recup = 
//             (std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
//           for(unsigned int j=0 ; j<(*recup).size() ; j++) {
//             int remoteID  = (*recup)[j].first;
//             if(remoteID == recvID) {
//               (*recup)[j].second = precv;
//             }
//           }  
//         }
//         AP_free(msg);
//       }		   
   
//     }
//     AP_check_sends(AP_WAITALL);
//     MPI_Barrier(MPI_COMM_WORLD);
//     delete [] sendcounts;
//   }
// #endif

  
  // -------------------------------------------------------------------
  /*! \brief Create edge correspondance and assure coherent orientation \ingroup parallel */ 
  void E_createInfoInterface(pMesh mesh, pMeshDataId tagEdge)
  {

#ifdef PARALLEL
    int nproc = 1;
#endif
    int myrank = 0;

    size_t nbConnections = 0;

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

    // --- Delete previous information ---
    
    EIter eit = M_edgeIter(mesh);
    while ( pEdge pe = EIter_next(eit) ) {
      void * tmp_ptr;
      if( EN_getDataPtr((pEntity) pe, tagEdge, &tmp_ptr) ) {
        EN_deleteData((pEntity) pe, tagEdge);
      }
    }
    EIter_delete(eit);


    // --- establish local connections and assure coherent orientation --- 
    // 
    //     * coherent reorientation on other partitions assumes coherent orientation here 
    //     * all periodic edges have same connections across partititions
    //       hence subsequent reorientation should remain coherent
    
    eit = M_edgeIter(mesh);
    
    while ( pEdge pe = EIter_next(eit) ) {

      std::vector<edge_comm> remote;
      
      int id0 = EN_id((pEntity) E_vertex(pe,0));
      int id1 = EN_id((pEntity) E_vertex(pe,1));

      std::pair<int,int> lNodeID(std::min(id0,id1),
                                 std::max(id0,id1));
      
      if (E_isPotentialInterface(mesh,pe,remote)) {
        
        std::vector<edge_comm >::const_iterator riter = remote.begin();


        pVertex p1Local = E_vertex(pe,0);
        pVertex p2Local = E_vertex(pe,1);
        

        for (;riter!=remote.end();++riter) {
        
          int distProc = riter->distProc;

          if (distProc == myrank) {
            
            edge_comm comm = *riter;
            
            pVertex p1 = comm.p1;
            pVertex p2 = comm.p2;
            
            pEdge remoteEdge = E_exist(p1,p2);
            
            if (remoteEdge) {
              

              // get points again, since alignment may be different from comm 

              pVertex p1Remote = E_vertex(remoteEdge,0);
              pVertex p2Remote = E_vertex(remoteEdge,1);

              // verify orientation 

              int orientation = 0;
              if (V_corresponds(p1Remote,p1Local) && V_corresponds(p2Remote,p2Local)) orientation =  1;
              if (V_corresponds(p1Remote,p2Local) && V_corresponds(p2Remote,p1Local)) orientation = -1;
                
              if (orientation == 0) {
                MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                            "Got a locally periodic edge connection from %d-%d to %d-%d which does not correspond nodewise",
                                            EN_id((pEntity) p1),EN_id((pEntity) p2),
                                            EN_id((pEntity) p1Local),EN_id((pEntity) p2Local));
              }
              
              std::pair<int,int> rNodeID(std::min(EN_id((pEntity) p1Remote),EN_id((pEntity) p2Remote)),
                                         std::max(EN_id((pEntity) p1Remote),EN_id((pEntity) p2Remote))); 
  
            
              if (rNodeID < lNodeID) { 
                lNodeID = rNodeID; // remember to what we aligned the edge 
                if (orientation == -1) {
                  int dir = E_align(pe,p2Local,p1Local);
                  if (dir != -1) std::cout << "misalignment of edges " << std::endl;
                }
              }
              
              void* tmp_ptr;
              std::multimap<int,pEdge>* list;
              
              if(EN_getDataPtr((pEntity) pe, tagEdge, (void**)&tmp_ptr)) {
                list = (std::multimap<int,pEdge>*) tmp_ptr;
              }
              else {
                list = new std::multimap<int,pEdge> ;
                EN_attachDataPtr((pEntity) pe,tagEdge,list);
              }

              list->insert(std::make_pair(myrank,remoteEdge));
              nbConnections++;
            }
          }
        }
      }
    }
    EIter_delete(eit);
    
#ifdef PARALLEL
    
    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    
    // --- send edge info ---

    eit = M_edgeIter(mesh);
    while ( pEdge pe = EIter_next(eit) ) {

      std::vector<edge_comm> remote;

      if (E_isPotentialInterface(mesh,pe,remote)) {

        std::vector<edge_comm >::const_iterator riter = remote.begin();

        for (;riter!=remote.end();++riter) {
        
          int distProc = riter->distProc;

          if (distProc != myrank) {
            void *buf = AP_alloc(distProc,446,sizeof(edge_comm));
            edge_comm *castbuf = (edge_comm *) buf;
            castbuf->p1  = riter->p1;
            castbuf->p2  = riter->p2;
            castbuf->id1 = riter->id1;
            castbuf->id2 = riter->id2;
            castbuf->pe = pe;
            AP_send(buf);
            sendcounts[distProc]++;
          }
        }
      }
    }
    EIter_delete(eit);
    
    // --- synchronise sends --- 
    
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
    //AP_flush();

    // --- Receive edge info ---
    
    int message=0,count;
    while (!AP_recv_count(&count) ||message<count) {

      void *msg;
      int from;
      int tag;
      int size;
      int rc;

      rc=AP_recv(MPI_ANY_SOURCE,446, AP_BLOCKING|AP_DROPOUT,
                 &msg, &size, &from, &tag);
      
      if (rc) {
        message++;
        edge_comm * comm = (edge_comm*) msg;
        pEdge pe = E_exist(comm->p1,comm->p2);
        if(pe) {
        
          int lowestRank = myrank;

          // avoid connections between two edges that span full periodicity
          // ie. periodic nodes n1 and n3 while both edges n1-n2 and n2-n3 exist

          int id1 = comm->id1;
          int id2 = comm->id2;

          if (E_corresponds(pe,id1,id2)) {
            
            std::multimap<int,pEdge>* list;
            if(EN_getDataPtr((pEntity) pe, tagEdge, (void**)&list)) {
              lowestRank = std::min(list->begin()->first,lowestRank);
            }
            else {
              list = new std::multimap<int,pEdge>;
              EN_attachDataPtr((pEntity) pe, tagEdge, (void*)list);
            }
            
            list->insert(std::pair<int,pEdge>(from,comm->pe));
            if (from < lowestRank) {
              int dir = E_align(pe,comm->p1,comm->p2);
            }
          }
        }
        AP_free(msg);
      }
    }

    // --- synchronise receives --- 
    
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // --- clean ship ---
    
    delete [] sendcounts;
    
#endif
    
    // --- verify !! --- 
    
    MDB_CommCheck cc;
    exchangeDataOnEdges(mesh,cc);

  }

  // -------------------------------------------------------------------
  /*! \brief Establish face to face correspondance and assure coherent orientation \ingroup parallel */ 
  //
  
  void F_createInfoInterface(pMesh mesh, pMeshDataId tagFace)
  {

#ifdef PARALLEL
    int nproc  = 1;
#endif
    int myrank = 0;
    size_t nbConnections = 0;
    
#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

    // --- delete previous information ---

    FIter fit = M_faceIter(mesh);
    while ( pFace pf = FIter_next(fit) ) {
      void * tmp_ptr;
      if( EN_getDataPtr((pEntity) pf, tagFace, &tmp_ptr) ) {
        EN_deleteData((pEntity) pf, tagFace);
      }
    }
    FIter_delete(fit);

    // --- first establish local communication and coherent orientation ---
    //     * coherent reorientation on other partitions assumes coherent orientation here 
    //     * all periodic faces have same connections across partititions
    //       hence subsequent reorientation should remain coherent
    
    
    fit = M_faceIter(mesh);
    while ( pFace pf = FIter_next(fit) ) {
      int distProc = -1;

      std::vector<pVertex> distVt;
      std::vector<face_comm> rd;
      
      if ( F_isPotentialInterface(mesh, pf, rd) ) {
        
        std::vector<face_comm>::iterator rIter = rd.begin();

        std::vector<int> lNodeID;
        for (int i=0;i<F_numVertices(pf);i++) lNodeID.push_back(EN_id((pEntity) F_vertex(pf,i)));
        
        
        for (;rIter!=rd.end();++rIter) {
          
          distProc = rIter->distProc;
          
          if (distProc == myrank) {
            
            face_comm& comm = *rIter;
            
            pVertex p1recv,p2recv,p3recv,p4recv;
            int nbNodes = (comm.p4 == NULL) ? 3:4;

            if (nbNodes != F_numVertices(pf)) {
              std::cout << "Non-coherent number of vertices " << nbNodes << " vs " << F_numVertices(pf) << std::endl;
            }
            
            p1recv = comm.p1;
            p2recv = comm.p2;
            p3recv = comm.p3;
            p4recv = comm.p4;
            
            assert(p1recv);assert(p2recv);assert(p3recv);
            // p4recv = nbNodes == 4 ? comm.p4 : NULL;
            
            pFace remoteFace = F_exist(p1recv,p2recv,p3recv,p4recv);
            
            if (remoteFace) {
            
              // compare ordered lists of nodes 
              
              std::vector<int> rNodeID;
              rNodeID.push_back(EN_id((pEntity) p1recv));
              rNodeID.push_back(EN_id((pEntity) p2recv));
              rNodeID.push_back(EN_id((pEntity) p3recv));
              if (nbNodes == 4) rNodeID.push_back(EN_id((pEntity) p4recv));
              
              // this is ok : only 1 periodic connection possible per face 

              if (rNodeID > lNodeID) {
                lNodeID = rNodeID; 
                F_align(remoteFace,p1recv,p2recv,p3recv,p4recv);
              }
              
              void* tmpptr;
              std::multimap<int,pFace>* list = NULL;
              
              if(EN_getDataPtr((pEntity) pf , tagFace,&tmpptr)){
                list = (std::multimap<int,pFace>*) (tmpptr);
              } else {
                list = new std::multimap<int,pFace>;
                EN_attachDataPtr((pEntity) pf, tagFace, list);
              }
              list->insert(std::make_pair(myrank,remoteFace));
              nbConnections++;
            }
          }
        }
      }
    }
    
#ifdef PARALLEL

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;

    // --- send faces ---     

    fit = M_faceIter(mesh);
    while ( pFace pf = FIter_next(fit) ) {
      int distProc = -1;
      std::vector<pVertex> distVt;
      
      std::vector<face_comm> rd;
      
      if ( F_isPotentialInterface(mesh, pf, rd) ) {
        
        std::vector<face_comm>::iterator rIter = rd.begin();
        
        for (;rIter!=rd.end();++rIter) {
          
          distProc = rIter->distProc;
          
          if (distProc != myrank) {
            
            void *buf = AP_alloc(distProc,444,sizeof(face_comm));
            face_comm *castbuf = (face_comm *) buf;
            memcpy(buf,&(*rIter),sizeof(face_comm));
            castbuf->distProc = myrank;
            
            AP_send(buf);
            sendcounts[distProc]++;
          }
        }
      }
    }
    FIter_delete(fit);

    // --- synchronise sends --- 

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
    //AP_flush();
  
    // --- Receive face info ---
    int message=0, count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int from;
      int tag;
      int size;
      int rc;
      rc=AP_recv(MPI_ANY_SOURCE,444, AP_BLOCKING|AP_DROPOUT,
                 &msg, &size, &from, &tag);
      if (rc) {
        message++;
        
        face_comm * castbuf = (face_comm*) msg;
        
        pVertex p1recv = castbuf -> p1;
        pVertex p2recv = castbuf -> p2;
        pVertex p3recv = castbuf -> p3;
        pVertex p4recv = castbuf -> p4;

        int id1recv = castbuf->id1;
        int id2recv = castbuf->id2;
        int id3recv = castbuf->id3;
        int id4recv = castbuf->id4;
      
        assert(p1recv);assert(p2recv);assert(p3recv);
        // p4recv = nbNodes == 4 ? castbuf -> p4 : NULL;
      
        int nprocrecv = castbuf->distProc;
        assert(nprocrecv==from);
      
        pFace pface = F_exist(p1recv,p2recv,p3recv,p4recv);
        
        if (pface) {

          if (F_corresponds(pface,id1recv,id2recv,id3recv,id4recv)) {

            // lowest ranking processor is master
            
            void* tmpptr;
            std::multimap<int,pFace>* list = NULL;
            
            if( EN_getDataPtr((pEntity) pface , tagFace,&tmpptr)){
              list = (std::multimap<int,pFace>*) (tmpptr);
              if (from < myrank && from < list->begin()->first) {
                F_align(pface,p1recv,p2recv,p3recv,p4recv);
              }
            } 
            else {
              if (from < myrank) F_align(pface,p1recv,p2recv,p3recv,p4recv);
              list = new std::multimap<int,pFace>;
              EN_attachDataPtr((pEntity) pface, tagFace, list);
            }
            list->insert(std::make_pair(nprocrecv,castbuf->pf));
          }
        }
        AP_free(msg);
      }
    }

    // --- synchronise receives --- 

    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);

    // --- clean ship --- 

    delete [] sendcounts;

#endif
    
    // --- verify !! --- 

    MDB_CommCheck cc;
    exchangeDataOnFaces(mesh,cc);
  }
  
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // #ifdef PARALLEL
  // void UpdateIDGlobal(pMesh mesh, int IdGlobal)
  // {
  //   int mysize,myrank;
  //   MPI_Comm_size(MPI_COMM_WORLD, &mysize);
  //   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //   // send pointers
  //   int *sendcounts = new int[mysize];
  //   for(int i=0; i<mysize; i++) sendcounts[i] = 0;
  
  //   pMeshDataId tagVertex = MD_lookupMeshDataId("RemotePoint");
  //   pMeshDataId tagGlob = MD_lookupMeshDataId("IdGlobal");

  //   VIter vit = M_vertexIter(mesh);
  //   pVertex pv;
  //   while ((pv = VIter_next(vit))) {
  //     void *temp_ptr; 
  //     int isInterface = EN_getDataPtr((pEntity) pv , tagVertex, &temp_ptr);
    
  //     if(isInterface) {
  //       std::vector<std::pair<int , MDB_Point*> > *recup = (std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
  //       int minproc = mysize + 10;
  //       for(unsigned int j=0 ; j<(*recup).size() ; j++) {
  //         minproc = std::min(minproc,(*recup)[j].first);
  //       }
  //       minproc = std::min(myrank,minproc);
  //       if(minproc!=myrank) continue;
  //       int numGlobal = EN_id((pEntity) pv) + IdGlobal + 1;
  //       EN_attachDataInt((pEntity) pv, tagGlob, numGlobal);
  //       for(unsigned int j=0 ; j<(*recup).size() ; j++) {
  //         int remoteID  = (*recup)[j].first;
  //         assert(remoteID != myrank);
  //         void *buf = AP_alloc(remoteID,444,sizeof(point_comm));
  //         point_comm *castbuf = (point_comm *) buf;
  //         castbuf->p         = (*recup)[j].second;
  //         castbuf->nID       =  myrank;
  //         castbuf->numGlobal =  numGlobal;
  //         AP_send(buf);
  //         sendcounts[remoteID]++;
  //       }
  //     }
  //   }
  //   VIter_delete(vit);

  //   AP_check_sends(AP_NOFLAGS);
  //   AP_reduce_nsends(sendcounts);
  //   AP_flush();

  //   // receive pointers
  //   int message=0;
  //   int count;
  //   while (!AP_recv_count(&count) || message<count) {
  //     void *msg;
  //     int  from;
  //     int  tag;
  //     int  size;
  //     int  recv;
  //     recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
  // 		   &msg, &size, &from, &tag);
  //     if(recv) {
  //       message++;
  //       point_comm * castbuf = (point_comm*) msg;
  //       pVertex precv (0);
  //       precv = castbuf -> p;
  //       assert(precv);
  //       //int recvID = castbuf -> nID;
  //       int numrecv = castbuf -> numGlobal;
      
  //       // DEBUG
  //       void *temp_ptr; 
  //       int isInterface = EN_getDataPtr((pEntity) precv, tagVertex, &temp_ptr);
  //       assert(isInterface);
      
  //       EN_attachDataInt((pEntity) precv, tagGlob, numrecv);        
  //       AP_free(msg);
  //     }
    
  //   }
  //   AP_check_sends(AP_WAITALL);
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   delete [] sendcounts;
  // }
  // #endif

  // -------------------------------------------------------------------

}
