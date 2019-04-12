// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Cecile Dobrzynski, Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#ifdef PARALLEL
#include "MeshDataBaseParallelInterface.h"
#endif
// #include "NullModel.h"
#include "MeshDataBaseComm.h"


#include "assert.h"

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#endif

// for memcpy  (gcc 4.3) 
#include <cstring>

namespace MAd {

  void exchangeDataOnVertices (pMesh m, MDB_DataExchanger &de )
  { 
  
    int mysize = 1;
    int myrank = 0;

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
#endif    
    
    int *sendcounts = new int[mysize];
    for(int i=0;i<mysize;i++)sendcounts[i]=0;
  
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    VIter vit = M_vertexIter(m);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      
      if(isInterface) {
        const std::vector<std::pair<int , pVertex> > *recup = 
          (const std::vector<std::pair<int , pVertex> > *) temp_ptr;
        
        for(unsigned int j=0 ; j<(*recup).size() ; j++) {    
          
          
          
          int sizebuf;      
          int iProc = (*recup)[j].first;
          pVertex remote = (*recup)[j].second;
    
          void *buf = de.sendData ((pEntity) pv, iProc, sizebuf );
          if (buf)
          {
            if (iProc == myrank)
            {
              de.receiveData ((pEntity) remote, myrank,buf);
              free(buf);
            }
            else
            {
#ifdef PARALLEL
              char *msg = (char*)AP_alloc(iProc,de.tag(),sizeof(pEntity*)+sizebuf);
              memcpy(msg,&remote,sizeof(pEntity*));
              memcpy(&msg[sizeof(pEntity*)],buf,sizebuf);
              free(buf);
              AP_send(msg);
              sendcounts[iProc] ++;
#else
              throw;
#endif
            }
          } 
        }
      }
    }      
    VIter_delete(vit);  
#ifdef PARALLEL 
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    int message=0;
    int count;
  
    while (!AP_recv_count(&count) || message<count) 
      {
        void *msg;
        int from;
        int tag;
        int size;
        int rc;
        rc=AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                   &msg, &size, &from, &tag);
        if (rc) 
          {
            message++;
            char * tmp = (char *) msg;
            pEntity * pv = (pEntity *) &tmp[0];
#ifdef DEBUG	  
            assert(*pv);
#endif
            de.receiveData (*pv,from, &tmp[sizeof(pEntity*)]);
            AP_free(msg);
          }
      }    

    AP_check_sends(AP_WAITALL); 
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    delete [] sendcounts;
  }

  // ----------------------------------------------------------------------------
  
  void exchangeDataOnEdges (pMesh m, MDB_DataExchanger &de )
  { 
    int mysize = 1;
    int myrank = 0;

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif    
  
    int *sendcounts = new int[mysize];
    for(int i=0;i<mysize;i++)sendcounts[i]=0;
    
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    EIter eit = M_edgeIter(m);
    pEdge pe;  
    while (pe = EIter_next(eit)) {

      void* tmpptr = NULL;
      if (EN_getDataPtr(pe,tagData,&tmpptr)) {
        
        std::multimap<int,pEdge>* connections = (std::multimap<int,pEdge>*) tmpptr;
        std::multimap<int,pEdge>::iterator cIter = connections->begin();
        
        
        for (;cIter!=connections->end();++cIter) {
          
          int iProc = cIter->first;
          pEdge remote = cIter->second;
          int sizebuf = 0;
          void *buf = de.sendData ((pEntity) pe, iProc, sizebuf );
          size_t sendSize = sizebuf + sizeof(pEntity*);
  
          if (iProc == myrank) de.receiveData ((pEntity) cIter->second,myrank,buf);
          else {
#ifdef PARALLEL
            char *msg = (char*)AP_alloc(iProc,de.tag(),sendSize);
            memcpy(&msg[0],&remote,sizeof(pEntity));
            memcpy(&msg[sizeof(pEntity)],buf,sizebuf);
            AP_send(msg);
            sendcounts[iProc] ++;
#else
            throw;
#endif
          }
          free(buf); 
        }
      }
    }

    EIter_delete(eit);
              

#ifdef PARALLEL 
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    int message=0;
    int count;
  
    while (!AP_recv_count(&count) || message<count) 
      {
        void *msg;
        int from;
        int tag;
        int size;
        int rc;
        rc=AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                   &msg, &size, &from, &tag);
        if (rc) 
          {
            message++;
            char * tmp = (char *) msg;
            pEdge pe;
            memcpy(&pe,&tmp[0],sizeof(pEntity));
            de.receiveData ((pEntity) pe,from, &tmp[sizeof(pEntity*)]);
            AP_free(msg);
          }
      }    
    AP_check_sends(AP_WAITALL); 
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete [] sendcounts;
  }

  // ----------------------------------------------------------------------------

  void exchangeDataOnFaces (pMesh m, MDB_DataExchanger &de )
  { 
    
    int mysize = 1;
    int myrank = 0;

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif    
  
    int *sendcounts = new int[mysize];
    for(int i=0;i<mysize;i++)sendcounts[i]=0;
    
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    FIter fit = M_faceIter(m);
    pFace pf;  
    while (pf = FIter_next(fit)) {

      void* tmpptr = NULL;
      if (EN_getDataPtr((pEntity) pf,tagData,&tmpptr)) {
        
        std::multimap<int,pFace>* connections = (std::multimap<int,pFace>*) tmpptr;
        std::multimap<int,pFace>::iterator cIter = connections->begin();
        
        
        for (;cIter!=connections->end();++cIter) {
          
          int iProc = cIter->first;
          pFace remote = cIter->second;
          int sizebuf = 0;
          void *buf = de.sendData ((pEntity) pf, iProc, sizebuf );
          size_t sendSize = sizebuf + sizeof(pEntity*);
  
          if (iProc == myrank) de.receiveData ((pEntity) cIter->second,myrank,buf);
          else {
            
#ifdef PARALLEL
            char *msg = (char*)AP_alloc(iProc,de.tag(),sendSize);
            memcpy(&msg[0],&remote,sizeof(pEntity));
            memcpy(&msg[sizeof(pEntity)],buf,sizebuf);
            AP_send(msg);
            sendcounts[iProc] ++;
#else
            throw;
#endif
          }
          free(buf); 
        }
      }
    }
    FIter_delete(fit);

#ifdef PARALLEL 
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    int message=0;
    int count;
  
    while (!AP_recv_count(&count) || message<count) 
      {
        void *msg;
        int from;
        int tag;
        int size;
        int rc;
        rc=AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                   &msg, &size, &from, &tag);
        if (rc) 
          {
            message++;
            char * tmp = (char *) msg;
            pEdge pe;
            memcpy(&pe,&tmp[0],sizeof(pEntity));
            de.receiveData ((pEntity) pe,from, &tmp[sizeof(pEntity*)]);
            AP_free(msg);
          }
      }    
    AP_check_sends(AP_WAITALL); 
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete [] sendcounts;
  }


//   void exchangeDataOnEdges_old (pMesh m, MDB_DataExchanger &de )
//   { 

//     int mysize = 1;
//     int myrank = 0;
// #ifdef PARALLEL
//     MPI_Comm_size(MPI_COMM_WORLD, &mysize);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
// #endif    
  
//     int *sendcounts = new int[mysize];
//     for(int i=0;i<mysize;i++)sendcounts[i]=0;
  
//     pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
//     EIter eit = M_edgeIter(m);
//     pEdge pe;  
//     while ((pe = EIter_next(eit))) {
// #ifdef PARALLEL
//       int distProc = -1;
//       std::vector<pVertex> distVt;

//       /// E_isInterface should be replaced by correct version or by direct use of the remotePoint
//         if ( E_isInterface(m,pe,&distProc,&distVt) ) {
     
//           //		void * tmp_ptr;
//           //    if( EN_getDataPtr((pEntity) pe, tagData, &tmp_ptr) ){
// #endif
//           pVertex p1 = pe->p1;
//           pVertex p2 = pe->p2;
//           assert(p1); assert(p2);
//           void *temp_ptr1,*temp_ptr2; 
//           EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
//           EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);       
//           const std::vector<std::pair<int , pVertex> > *recup1 = 
//             (const std::vector<std::pair<int , pVertex> > *) temp_ptr1;
//           const std::vector<std::pair<int , pVertex> > *recup2 = 
//             (const std::vector<std::pair<int , pVertex> > *) temp_ptr2;
//           int size = (*recup1).size();
//           assert(size);
//           int tab[1024];
//           for(int i=0 ; i< size ; i++) tab[i] = (*recup1)[i].first;
	
//           for(unsigned int j=0 ; j<(*recup2).size() ; j++) {    
//             int iProc = (*recup2)[j].first;
//             int i;
//             for(i=0 ; i<size ; i++) {
//               if(iProc == tab[i]) break;
//             }
//             if(i < size) {
//               assert(iProc == tab[i]);
	    
//               int sizebuf;
//               pVertex remote1 = (*recup1)[i].second;
//               pVertex remote2 = (*recup2)[j].second;
//               void *buf = de.sendData ((pEntity) pe, iProc, sizebuf );
//               if (buf)
//                 {
//                   if (iProc == myrank)
//                     {
//                       puts("*****************************pas parallele");
//                       char *msg = (char*)malloc(2*sizeof(pEntity)+sizebuf);
//                       memcpy(msg,&remote1,sizeof(pEntity));
//                       memcpy(&msg[sizeof(pEntity)],&remote2,sizeof(pEntity));
//                       memcpy(&msg[2*sizeof(pEntity)],buf,sizebuf);
//                       free(buf);
//                       MDB_VectorE ve = remote1->edges;
//                       MDB_VectorE::iterator it  = ve.begin();
//                       MDB_VectorE::iterator ite = ve.end();
//                       while(it!=ite){
//                         if((*it)->p1 == remote1 && (*it)->p2 == remote2)break;
//                         else if ((*it)->p2 == remote1 && (*it)->p1 == remote2) break;
//                         ++it;
//                       }
//                       assert(it!=ite);
//                       pEdge pe = *it;
//                       de.receiveData ((pEntity) pe, 
//                                       myrank,
//                                       msg);
//                       free(msg);
//                     }
//                   else
//                     {
// #ifdef PARALLEL
//                       char *msg = (char*)AP_alloc(iProc, 
//                                                   de.tag() , 
//                                                   2*sizeof(pEntity*)+sizebuf);			   
//                       memcpy(&msg[0],&remote1,sizeof(pEntity));
//                       memcpy(&msg[sizeof(pEntity)],&remote2,sizeof(pEntity));
//                       memcpy(&msg[2*sizeof(pEntity)],buf,sizebuf);
//                       free(buf);
//                       AP_send(msg);
//                       sendcounts[iProc] ++;
// #else
//                       throw;
// #endif
//                     } 
//                 } 
//             }
//           }
// #ifdef PARALLEL
//         }
// #endif
//     }      
//     EIter_delete(eit);  

// #ifdef PARALLEL 
//     AP_check_sends(AP_NOFLAGS);
//     AP_reduce_nsends(sendcounts);
  
//     int message=0;
//     int count;
  
//     while (!AP_recv_count(&count) || message<count) 
//       {
//         void *msg;
//         int from;
//         int tag;
//         int size;
//         int rc;
//         rc=AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
//                    &msg, &size, &from, &tag);
//         if (rc) 
//           {
//             message++;
//             char * tmp = (char *) msg;
//             pVertex v1;
//             pVertex v2;
//             memcpy(&v1,&tmp[0],sizeof(pEntity));
//             memcpy(&v2,&tmp[sizeof(pEntity)],sizeof(pEntity));

//             MDB_VectorE &ve = v1->edges;
//             MDB_VectorE::iterator it  = ve.begin();
//             MDB_VectorE::iterator ite = ve.end();
//             while(it!=ite){
//               if((*it)->p2 == (v2) || (*it)->p1 == (v2)) 
//                 {
//                   assert((*it)->p2 == (v1) || (*it)->p1 == (v1));
//                   break;
//                 }    
//               ++it;
//             }
//             // MISTAKE HERE, CECILE, YOU SHOULD VERIFY THAT !!
//             if(it==ite)
//               {
//                 // 	    printf("the edge does not exist here...\n");
//                 //	    throw;
//               }
//             else
//               {
//                 pEdge pe = *it;
	      
//                 /*
//                   Little test : if edges are reversed on the
//                   different sides of the interface then reverse
//                   the ones on the smallest processor
//                 */		
//                 if (v1 == pe->p2 && from < myrank)
//                   {
//                     pe->p1 = v1;
//                     pe->p2 = v2;
//                   } 
//                 // ------- END OF TEST ------------------------


//                 de.receiveData ((pEntity) pe,from, &tmp[2*sizeof(pEntity*)]);
//               }
//             AP_free(msg);
//           }
//       }    
//     AP_check_sends(AP_WAITALL); 
//     MPI_Barrier(MPI_COMM_WORLD);
// #endif
//     delete [] sendcounts;
//   }


  // ---------------------------------------------------------------------
  /*
    The parallel faces are supposed to be classified on a geometric entity 
    of dimension 2 with the same tag as the unique volume entity which 
    supports the face on this partition. This is done when a mesh is 
    partitioned with japp.
  */
//   void exchangeDataOnFaces_old (pMesh m, MDB_DataExchanger &de )
//   {
//     exchangeDataOnTriangles(m, de);
//     exchangeDataOnQuads(m, de);
//   }

//   void exchangeDataOnTriangles (pMesh m, MDB_DataExchanger &de )
//   { 
//     int mysize = 1;
//     int myrank = 0;
    
// #ifdef PARALLEL
//     MPI_Comm_size(MPI_COMM_WORLD, &mysize);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
// #endif    
  
//     // --------------------------------------------------------
//     // --- SEND -----------------------------------------------
//     // --------------------------------------------------------
  
//     int *sendcounts = new int[mysize];
//     for(int i=0;i<mysize;i++) sendcounts[i]=0;
  
//     pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
//     //  FIter fit = M_faceIter(m);
//     MDB_FIter fit(&m->triangles);
//     pFace pface;  
//     //  while ((pface = FIter_next(fit))) {
//     while (( pface = fit.next() )) {
//       pVertex nod[3];
//       pface->getNodes(nod);
//       pVertex p1 = nod[0];
//       pVertex p2 = nod[1];
//       pVertex p3 = nod[2];
//       void *temp_ptr1,*temp_ptr2,*temp_ptr3; 
//       int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
//       int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
//       int isInterface3 = EN_getDataPtr((pEntity) p3 , tagData, &temp_ptr3);
//       if(isInterface1 && isInterface2 && isInterface3) {
//         const std::vector<std::pair<int , pVertex> > *recup1 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr1;
//         const std::vector<std::pair<int , pVertex> > *recup2 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr2;
//         const std::vector<std::pair<int , pVertex> > *recup3 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr3;
//         int size = (*recup1).size();
//         int *tab=new int[size];
//         for(int i=0 ; i< size ; i++) tab[i] = (*recup1)[i].first;
      
//         int size2 = (*recup2).size();
//         int *tab2=new int[size2];
//         for(int i=0 ; i< size2 ; i++) tab2[i] = (*recup2)[i].first;
      
//         for(unsigned int k=0 ; k<(*recup3).size() ; k++) {  
//           // the 3 vertices could belong to 3 identical procs at the same time and the face 
//           // belongs to only 2 of those. In this case 2 packets are sent. We have to 
//           // ignore the 'wrong' packet at reception. (Done)
//           int iProc = (*recup3)[k].first;
//           int i;
//           for(i=0 ; i<size ; i++) {
//             if(iProc == tab[i]) break;
//           }
//           int j;
//           for(j=0 ; j<size2 ; j++) {
//             if(iProc == tab2[j]) break;
//           }
//           if(i < size && j < size2) {
//             assert(tab[i]==iProc);
//             assert(tab2[j]==iProc);
	      
//             int sizebuf;
          
//             pVertex remote1 = (*recup1)[i].second;
//             pVertex remote2 = (*recup2)[j].second;
//             pVertex remote3 = (*recup3)[k].second;
	  
//             void *buf = de.sendData ((pEntity) pface, iProc, sizebuf );
//             if (buf) {
//               if (iProc == myrank) {
//                 puts("*****************************pas parallele");
//                 char *msg = (char*)malloc(3*sizeof(pEntity*)+sizebuf);
//                 memcpy(msg,&remote1,sizeof(pEntity*));
//                 memcpy(&msg[sizeof(pEntity*)],&remote2,sizeof(pEntity*));
//                 memcpy(&msg[2*sizeof(pEntity*)],&remote3,sizeof(pEntity*));
//                 memcpy(&msg[3*sizeof(pEntity*)],buf,sizebuf);
//                 free(buf);
//                 MDB_ListF vf;
//                 remote1->getTriangles(vf);
//                 MDB_ListF::iterator it  = vf.begin();
//                 MDB_ListF::iterator ite = vf.end();
//                 while(it!=ite){
//                   pVertex nod[3];
//                   (*it)->getNodes(nod);
//                   if((nod[0] == remote1 || nod[1] == remote1 || nod[2] == remote1)&& 
//                      (nod[0] == remote2 || nod[1] == remote2 || nod[2] == remote2)&&
//                      (nod[0] == remote3 || nod[1] == remote3 || nod[2] == remote3))break;
//                   ++it;
//                 }
//                 assert(it!=ite);
//                 pFace pface = *it;
//                 de.receiveData ((pEntity) pface, 
//                                 myrank,
//                                 msg);
//                 free(msg); 
//               }
//               else {
// #ifdef PARALLEL
//                 char *msg = (char*)AP_alloc(iProc, 
//                                             de.tag() , 
//                                             3*sizeof(pEntity*)+sizebuf);
//                 memcpy(&msg[0],&remote1,sizeof(pEntity*));
//                 memcpy(&msg[sizeof(pEntity*)],&remote2,sizeof(pEntity*));
//                 memcpy(&msg[2*sizeof(pEntity*)],&remote3,sizeof(pEntity*));
//                 memcpy(&msg[3*sizeof(pEntity*)],buf,sizebuf);
//                 free(buf);
//                 AP_send(msg);
//                 sendcounts[iProc]++;
// #else
//                 throw;
// #endif
//               } 
//             }
//           }
//         }
//         delete []tab;
//         delete []tab2;
//       }
//     }      
//     //  FIter_delete(fit);  
// #ifdef PARALLEL
//     // #ifdef MDB_DEBUG
//     //   MPI_Barrier( MPI_COMM_WORLD );
//     //   AP_check_sends(AP_WAITDEFER);
//     // #else
//     AP_check_sends(AP_NOFLAGS);
//     // #endif
//     AP_reduce_nsends(sendcounts);
  
//     // --------------------------------------------------------
//     // --- RECEIVE --------------------------------------------
//     // --------------------------------------------------------
  
//     int message=0;
//     int count;

//     while (!AP_recv_count(&count) || message<count) 
//       {
//         void *msg;
//         int from;
//         int tag;
//         int size;
//         int rc;
//         rc=AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
//                    &msg, &size, &from, &tag);
//         if (rc) 
//           {
//             message++;
//             char * tmp = (char *) msg;
//             pVertex * pv1 = (pVertex *) &tmp[0];
//             pVertex * pv2 = (pVertex *) &tmp[sizeof(pEntity*)];
//             pVertex * pv3 = (pVertex *) &tmp[2*sizeof(pEntity*)];
//             MDB_ListF vf;
//             (*pv1)->getTriangles(vf);
//             MDB_ListF::iterator it  = vf.begin();
//             MDB_ListF::iterator ite = vf.end();
//             while(it!=ite){
//               pVertex nod[3];
//               (*it)->getNodes(nod);
//               if((nod[0] == (*pv2) || nod[1] == (*pv2) || nod[2] == (*pv2))&&
//                  (nod[0] == (*pv3) || nod[1] == (*pv3) || nod[2] == (*pv3))) {
//                 assert((nod[0] == (*pv1) || nod[1] == (*pv1) || nod[2] == (*pv1)));
//                 break;
//               }
//               ++it;
//             }
//             // this can happen: see the comment above in the send part.
//             if ( it == ite ){
//               AP_free(msg);
//             }
//             else {
//               pFace pface = *it;
//               de.receiveData ((pEntity) pface,from, &tmp[3*sizeof(pEntity*)]);
//               AP_free(msg);
//             }
//           }
//       }    
//     AP_check_sends(AP_WAITALL); 
//     MPI_Barrier(MPI_COMM_WORLD);
// #endif
//     delete [] sendcounts;
//   }

//   void exchangeDataOnQuads (pMesh m, MDB_DataExchanger &de )
//   { 
//     int mysize = 1;
//     int myrank = 0;
    
// #ifdef PARALLEL
//     MPI_Comm_size(MPI_COMM_WORLD, &mysize);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
 
// #endif    
  
//     // --------------------------------------------------------
//     // --- SEND -----------------------------------------------
//     // --------------------------------------------------------
  
//     int *sendcounts = new int[mysize];
//     for(int i=0;i<mysize;i++) sendcounts[i]=0;
  
//     pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
//     //  FIter fit = M_faceIter(m);
//     MDB_QIter fit(&m->quads);
//     pFace pface;  
//     //  while ((pface = FIter_next(fit))) {
//     while ((pface = fit.next())) {
//       pVertex nod[4];
//       pface->getNodes(nod);
//       pVertex p1 = nod[0];
//       pVertex p2 = nod[1];
//       pVertex p3 = nod[2];
//       pVertex p4 = nod[3];
//       void *temp_ptr1,*temp_ptr2,*temp_ptr3,*temp_ptr4; 
//       int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
//       int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
//       int isInterface3 = EN_getDataPtr((pEntity) p3 , tagData, &temp_ptr3);
//       int isInterface4 = EN_getDataPtr((pEntity) p4 , tagData, &temp_ptr4);
//       if(isInterface1 && isInterface2 && isInterface3 && isInterface4) {
//         const std::vector<std::pair<int , pVertex> > *recup1 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr1;
//         const std::vector<std::pair<int , pVertex> > *recup2 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr2;
//         const std::vector<std::pair<int , pVertex> > *recup3 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr3;
//         const std::vector<std::pair<int , pVertex> > *recup4 = 
//           (const std::vector<std::pair<int , pVertex> > *) temp_ptr4;
//         int size = (*recup1).size();
//         int *tab=new int[size];
//         for(int i=0 ; i< size ; i++) tab[i] = (*recup1)[i].first;
      
//         int size2 = (*recup2).size();
//         int *tab2=new int[size2];
//         for(int i=0 ; i< size2 ; i++) tab2[i] = (*recup2)[i].first;

//         int size3 = (*recup3).size();
//         int *tab3=new int[size3];
//         for(int i=0 ; i< size3 ; i++) tab3[i] = (*recup3)[i].first;
      
//         for(int l=0 ; l<(int)((*recup4).size()) ; l++) {  
//           // this test is not robust: the 4 vertices could belong to 4 identical procs at the same
//           // time and the face belongs to only 2 of those. In this case 2 packets are sent. We have to 
//           // ignore the 'wrong' packet at reception. (Done)
//           int iProc = (*recup4)[l].first;
//           int i;
//           for(i=0 ; i<size ; i++) {
//             if(iProc == tab[i]) break;
//           }
//           int j;
//           for(j=0 ; j<size2 ; j++) {
//             if(iProc == tab2[j]) break;
//           } 
//           int k;
//           for(k=0 ; k<size3 ; k++) {
//             if(iProc == tab3[k]) break;
//           } 
//           if(i < size && j < size2 && k < size3) {
//             assert(tab[i]==iProc);
//             assert(tab2[j]==iProc);
//             assert(tab3[k]==iProc);

//             int dimGeomFace = GEN_type(pface->g);
//             if (dimGeomFace == 2) {
//               int tagFace = GEN_tag(pface->g); //std::cout <<"tagFace: "<<tagFace<<std::endl;
//               if (F_numRegions(pface)!=1) continue; //SL assert(F_numRegions(pface)==1);
//               int tagVol = GEN_tag(F_region(pface,0)->g); //std::cout <<"tagVol: "<<tagVol<<std::endl;
//               if ( tagFace == tagVol ) {
        
//                 int sizebuf;
        
//                 pVertex remote1 = (*recup1)[i].second;
//                 pVertex remote2 = (*recup2)[j].second;
//                 pVertex remote3 = (*recup3)[k].second;
//                 pVertex remote4 = (*recup4)[l].second;
        
//                 void *buf = de.sendData ((pEntity) pface, iProc, sizebuf );
//                 if (buf)
//                   {
//                     if (iProc == myrank)
//                       {
//                         puts("*****************************pas parallele");
//                         char *msg = (char*)malloc(4*sizeof(pEntity*)+sizebuf);
//                         memcpy(msg,&remote1,sizeof(pEntity*));
//                         memcpy(&msg[sizeof(pEntity*)],&remote2,sizeof(pEntity*));
//                         memcpy(&msg[2*sizeof(pEntity*)],&remote3,sizeof(pEntity*));
//                         memcpy(&msg[3*sizeof(pEntity*)],&remote4,sizeof(pEntity*));
//                         memcpy(&msg[4*sizeof(pEntity*)],buf,sizebuf);
//                         free(buf);
//                         MDB_ListF vf;
//                         remote1->getTriangles(vf);
//                         MDB_ListF::iterator it  = vf.begin();
//                         MDB_ListF::iterator ite = vf.end();
//                         while(it!=ite){
//                           pVertex nod[4];
//                           (*it)->getNodes(nod);
//                           if((nod[0] == remote1 || nod[1] == remote1 || nod[2] == remote1 || nod[3] == remote1)&& 
//                              (nod[0] == remote2 || nod[1] == remote2 || nod[2] == remote2 || nod[3] == remote2)&&
//                              (nod[0] == remote3 || nod[1] == remote3 || nod[2] == remote3 || nod[3] == remote3)&&
//                              (nod[0] == remote4 || nod[1] == remote4 || nod[2] == remote4 || nod[3] == remote4))break;
//                           ++it;
//                         }
//                         assert(it!=ite);
//                         pFace pface = *it;
//                         de.receiveData ((pEntity) pface, 
//                                         myrank,
//                                         msg);
//                         free(msg); 
//                       }
//                     else
//                       {
// #ifdef PARALLEL
//                         char *msg = (char*)AP_alloc(iProc, 
//                                                     de.tag() , 
//                                                     4*sizeof(pEntity*)+sizebuf);
//                         memcpy(&msg[0],&remote1,sizeof(pEntity*));
//                         memcpy(&msg[sizeof(pEntity*)],&remote2,sizeof(pEntity*));
//                         memcpy(&msg[2*sizeof(pEntity*)],&remote3,sizeof(pEntity*));
//                         memcpy(&msg[3*sizeof(pEntity*)],&remote4,sizeof(pEntity*));
//                         memcpy(&msg[4*sizeof(pEntity*)],buf,sizebuf);
//                         free(buf);
//                         AP_send(msg);
//                         sendcounts[iProc]++;
// #else
//                         throw;
// #endif
//                       } 
//                   }
//               }
//             } 
//           }
//         }
//         delete []tab;
//         delete []tab2;
//         delete []tab3;
//       }
//     }      
//     //  FIter_delete(fit);  
// #ifdef PARALLEL
//     // #ifdef MDB_DEBUG
//     //   //MPI_Waitall();
//     //   MPI_Barrier( MPI_COMM_WORLD );
//     //   AP_check_sends(AP_WAITALL);
//     // #else
//     AP_check_sends(AP_NOFLAGS);
//     // #endif
//     AP_reduce_nsends(sendcounts);
  
//     // --------------------------------------------------------
//     // --- RECEIVE --------------------------------------------
//     // --------------------------------------------------------
  
//     int message=0;
//     int count;
  
//     while (!AP_recv_count(&count) || message<count) 
//       {
//         void *msg;
//         int from;
//         int tag;
//         int size;
//         int rc;
//         rc=AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
//                    &msg, &size, &from, &tag);
//         if (rc) 
//           {
//             message++;
//             char * tmp = (char *) msg;
//             pVertex * pv1 = (pVertex *) &tmp[0];
//             pVertex * pv2 = (pVertex *) &tmp[sizeof(pEntity*)];
//             pVertex * pv3 = (pVertex *) &tmp[2*sizeof(pEntity*)];
//             pVertex * pv4 = (pVertex *) &tmp[3*sizeof(pEntity*)];
//             MDB_ListF vf;
//             (*pv1)->getTriangles(vf);
//             MDB_ListF::iterator it  = vf.begin();
//             MDB_ListF::iterator ite = vf.end();
//             while(it!=ite){
//               pVertex nod[4];
//               (*it)->getNodes(nod);
//               if((nod[0] == (*pv2) || nod[1] == (*pv2) || nod[2] == (*pv2) || nod[3] == (*pv2))&&
//                  (nod[0] == (*pv3) || nod[1] == (*pv3) || nod[2] == (*pv3) || nod[3] == (*pv3))&&
//                  (nod[0] == (*pv4) || nod[1] == (*pv4) || nod[2] == (*pv4) || nod[3] == (*pv4))) {
//                 assert((nod[0] == (*pv1) || nod[1] == (*pv1) || nod[2] == (*pv1) || nod[3] == (*pv1) ));
//                 break;
//               }
//               ++it;
//             }
//             // this can happen: see the comment above in the send part.
//             if ( it == ite ){
//               AP_free(msg);
//             }
//             else {
//               pFace pface = *it;
//               de.receiveData ((pEntity) pface,from, &tmp[4*sizeof(pEntity*)]);
//               AP_free(msg);
//             }
//           }
//       }    
//     AP_check_sends(AP_WAITALL); 
//     MPI_Barrier(MPI_COMM_WORLD);
// #endif
//     delete [] sendcounts;
//   }

} // End of namespace MAd
