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

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"

#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseParallelInterface.h"
#include "MeshDataBaseLoadBalancing.h"
#include "MeshDataBaseComm.h"

#include "assert.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>

namespace MAd {

  typedef std::pair<int, pVertex>              CommonpVertex_type;
  typedef std::vector< CommonpVertex_type  >  VectorOfCommonpVertex_type;
  // -------------------------------------------------------------------
  struct points_comm
  {
    pVertex p;     //pointeur in dest
    int nID;       //myrank
    int newproc;
  };

  // -------------------------------------------------------------------
  struct oldinterfaces_comm
  {
    int     oldnum;
    pVertex oldpv;                
  };
  // -------------------------------------------------------------------
  struct coor_comm
  {
    double   X,Y,Z;
    int tag,dim;
    int nproc;
    pVertex  psend;
  };  
  // -------------------------------------------------------------------
  struct newinterfaces_comm
  {
    int newproc;
  };  
  // -------------------------------------------------------------------
  struct points_comm_new
  {
    pVertex pdest;       //pointeur in dest
    int     destID;
    pVertex psend;
    int     sendID;       
  };

  // -------------------------------------------------------------------
  void ExchangeUpdatedInterfaceNodalInfo(pMesh mesh)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
 
    VIter vit = M_vertexIter(mesh);
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      void *temp_ptr;
      int isInterface = EN_getDataPtr((pEntity) pv, tagData, &temp_ptr);
    
      if(!isInterface) continue;
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      int tmp;
      int isChanged = EN_getDataInt((pEntity) pv, tagChange, &tmp);
      assert(isChanged);
      void *temp_ptr2;
      int is = EN_getDataPtr((pEntity) pv, tagRemote, &temp_ptr2);
      assert(is);
      VectorOfCommonpVertex_type *recup2 = 
        (VectorOfCommonpVertex_type  *) temp_ptr2;

      for(int j=0; j<(int)((*recup).size()); j++) {
        for( int i=0; i<(int)((*recup2).size()); i++) {
          int remoteID  = (*recup)[j].first;
          assert(remoteID != myrank);
          void *buf = AP_alloc(remoteID,444,sizeof(points_comm));
          points_comm *castbuf = (points_comm *) buf;
          castbuf->p     = (*recup)[j].second;
          castbuf->nID     = myrank;
          assert((*recup2)[i].first < nproc);
          castbuf->newproc   = (*recup2)[i].first;
          AP_send(buf);
          sendcounts[remoteID]++;
        }
      }
    }
    VIter_delete(vit);  
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        points_comm * castbuf = (points_comm*) msg;
        pVertex precv (NULL);
        precv = castbuf -> p;
        assert(precv);
        int nprocrecv = castbuf->newproc;
        assert(nprocrecv<nproc);
        int tmp; 
        int isChanged = EN_getDataInt((pEntity) precv, tagChange, &tmp);
        assert(isChanged);
        void *temp_ptr; 
        int is = EN_getDataPtr((pEntity) precv, tagRemote, &temp_ptr);
        assert(is);
        VectorOfCommonpVertex_type *recup = 
          (VectorOfCommonpVertex_type  *) temp_ptr;
        int * tabproc = new int[nproc];
        for(int i=0 ; i<nproc ; i++){
          tabproc[i] = 0;
        }
        for(unsigned int i=0 ; i<(*recup).size() ; i++){
          tabproc[(*recup)[i].first] = 1;
        }
        if(!tabproc[nprocrecv]){
          (*recup).push_back( CommonpVertex_type (nprocrecv,NULL));
        }
                
        AP_free(msg);
        delete []tabproc;
      }
    }
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;

  }

  // -------------------------------------------------------------------
  void MarkEltVertex(pMesh mesh, pMeshDataId tagElt)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Integer set to one if the vertex will change
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");

    // Data attached to each parallel vertex: distant proc 
    // together with distant pointer for the same vertex
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");

    // Like tagData but after the interface move
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");

  
    // --- Mark vertices ---
    int Dim = (mesh->tets.empty()) ? 2 : 3;
    if(Dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        int dest; 
        int migre = EN_getDataInt((pEntity) pface, tagElt, &dest);
        if(!migre) continue;    
    
        pVertex nod[3];
        nod[0] = F_vertex(pface,0);
        nod[1] = F_vertex(pface,1);
        nod[2] = F_vertex(pface,2);
        for(int i=0 ; i<3 ; i++) {
          EN_attachDataInt((pEntity) nod[i], tagChange, 1);
        }
      }
      FIter_delete(fit);
    } else {
      if(Dim!=3) throw;
      RIter rit = M_regionIter(mesh);
      pRegion pr;  
      while ((pr = RIter_next(rit))) {
        int dest; 
        int migre = EN_getDataInt((pEntity) pr, tagElt, &dest);
        if(!migre) continue;    
    
        pVertex nod[4];
        nod[0] = R_vertex(pr,0);
        nod[1] = R_vertex(pr,1);
        nod[2] = R_vertex(pr,2);
        nod[3] = R_vertex(pr,3);
        for(int i=0 ; i<4 ; i++) {
          EN_attachDataInt((pEntity) nod[i], tagChange, 1);
        }
      }
      RIter_delete(rit);
    }
  
    // --- Attach info about every future distant node to the marked nodes ---
    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    vit = M_vertexIter(mesh);
    while ((pv = VIter_next(vit))) {
      int tmp;
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv, tagData, &temp_ptr);
      if(!(isInterface || isChanged)) continue;
      if(!isChanged) EN_attachDataInt((pEntity) pv, tagChange, 1);
      /*vertex ball*/
      int ncompt = 0;
      int *tabproc=new int[nproc];
      for(int i=0 ; i<nproc ; i++){
        tabproc[i] = 0;
      }
      if(Dim==2) { 
        void *iter = 0;
        pPList ball = V_faces (pv);
        pFace pface;
        while( (pface =(pFace) PList_next(ball,&iter)) ) {
          int dest;
          int migre = EN_getDataInt((pEntity) pface, tagElt, &dest);
          if(migre) {
            tabproc[dest-1] = 1;
          } else {
            tabproc[myrank] = 1;
          }
        }
      } else {
        void *iter = 0;
        pPList ball = V_regions (pv);
        pRegion pr;
        while( (pr =(pRegion) PList_next(ball,&iter)) ) {
          int dest;
          int migre = EN_getDataInt((pEntity) pr , tagElt,&dest);
          if(migre) {
            tabproc[dest-1] = 1;
          } else {
            tabproc[myrank] = 1;
          }
        }
    
      }
      for(int i=0 ; i<nproc ; i++){
        ncompt += tabproc[i];
      }
      if(ncompt > 1) {
        VectorOfCommonpVertex_type  remoteUpdate;
        for(int i=0 ; i<nproc ; i++) {
          if(tabproc[i]) {
            (remoteUpdate).push_back(CommonpVertex_type(i,NULL));
          }
          EN_attachDataPtr((pEntity) pv , tagRemote, 
                           new VectorOfCommonpVertex_type (remoteUpdate));
        } 
      } else {
        int i=0;
        for(i=0 ; i<nproc ; i++) {
          if(tabproc[i] && (i != myrank)) break;  
        }  
        if(i==nproc) {
          assert(tabproc[myrank]);
          assert(isInterface);
          i = myrank;
        }
        VectorOfCommonpVertex_type  remoteUpdate;
        (remoteUpdate).push_back(CommonpVertex_type(i,NULL));
        EN_attachDataPtr((pEntity) pv , tagRemote, 
                         new VectorOfCommonpVertex_type (remoteUpdate));
      }
      delete []tabproc;
    }
    VIter_delete(vit);

    // --- Send and receive info at nodes for updated interface ---
    ExchangeUpdatedInterfaceNodalInfo(mesh);

  }

#ifdef DEBUG
  // -------------------------------------------------------------------
  void checkRemotePointer(pMesh mesh, pMeshDataId tagData ) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
 
    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      if(!isInterface) continue;
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      int sizeinterface = (*recup).size();

      for(int i=0 ; i<sizeinterface ; i++) {
        int remoteID = (*recup)[i].first;
        assert(remoteID!=myrank);
        void *buf = AP_alloc(remoteID,444,sizeof(oldinterfaces_comm));
        oldinterfaces_comm *castbuf = (oldinterfaces_comm *) buf; 
        castbuf->oldnum = myrank;
        assert((*recup)[i].second);
        castbuf->oldpv  = (*recup)[i].second;
        AP_send(buf);
        sendcounts[remoteID]++;       
      } 
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        oldinterfaces_comm * castbuf = (oldinterfaces_comm *) msg;
        pVertex pv    = castbuf->oldpv;
        int recvID    = castbuf->oldnum;
        void *temp_ptr; 
        int isInte = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
        assert(isInte);
        const VectorOfCommonpVertex_type *recup = 
          (const VectorOfCommonpVertex_type *) temp_ptr;
        unsigned int j=0;
        for( j=0 ; j<(*recup).size() ; j++) {
          int remoteID  = (*recup)[j].first;
          assert((*recup)[j].second);
          if(remoteID == recvID) break;
        }
        assert(j!= (*recup).size());  
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;



  }

  // -------------------------------------------------------------------
  void checkRemotePointer2(pMesh mesh, pMeshDataId tagData ) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
 
    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      if(!isInterface) continue;
      int tmp;
      int isChange=EN_getDataInt((pEntity) pv , tagChange, &tmp);
      if(!(isChange && tmp==1)) continue;
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      int sizeinterface = (*recup).size();

      for(int i=0 ; i<sizeinterface ; i++) {
        int remoteID = (*recup)[i].first;
        assert(remoteID!=myrank);
        void *buf = AP_alloc(remoteID,444,sizeof(oldinterfaces_comm));
        oldinterfaces_comm *castbuf = (oldinterfaces_comm *) buf; 
        castbuf->oldnum = myrank;
        castbuf->oldpv  = (*recup)[i].second;
        AP_send(buf);
        sendcounts[remoteID]++;       
      } 
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        oldinterfaces_comm * castbuf = (oldinterfaces_comm *) msg;
        pVertex pv    = castbuf->oldpv;
        int recvID    = castbuf->oldnum;
        void *temp_ptr; 
        int isInte = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
        assert(isInte);
        const VectorOfCommonpVertex_type *recup = 
          (const VectorOfCommonpVertex_type *) temp_ptr;
        unsigned int j=0;
        for( j=0 ; j<(*recup).size() ; j++) {
          int remoteID  = (*recup)[j].first;
          if(remoteID == recvID) break;
        }
        assert(j!= (*recup).size());  
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;



  }

  // -------------------------------------------------------------------
  void checkRemotePointerChange(pMesh mesh, pMeshDataId tagData,
                                pMeshDataId tagNew, pMeshDataId tagChange ) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
 
    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp;
      int isChange = EN_getDataInt((pEntity) pv , tagChange, &tmp);
      if(isChange) {
        void *temp_ptr; 
        int isInterface = EN_getDataPtr((pEntity) pv , tagNew, &temp_ptr);
        assert(isInterface);
        if(!isInterface) continue;
        VectorOfCommonpVertex_type *recup = 
          (VectorOfCommonpVertex_type *) temp_ptr;
        int sizeinterface = (*recup).size();
        for(int i=0 ; i<sizeinterface ; i++) {
          int remoteID = (*recup)[i].first;
          if(remoteID==myrank){
            assert(sizeinterface==1);
            continue;
          }
          void *buf = AP_alloc(remoteID,444,sizeof(oldinterfaces_comm));
          oldinterfaces_comm *castbuf = (oldinterfaces_comm *) buf; 
          castbuf->oldnum = myrank;
          castbuf->oldpv  = (*recup)[i].second;
          AP_send(buf);
          sendcounts[remoteID]++;       
        } 
    
      } else {
        void *temp_ptr; 
        int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
        if(!isInterface) continue;
        const VectorOfCommonpVertex_type *recup = 
          (const VectorOfCommonpVertex_type *) temp_ptr;
        int sizeinterface = (*recup).size();

        for(int i=0 ; i<sizeinterface ; i++) {
          int remoteID = (*recup)[i].first;
          assert(remoteID!=myrank);
          void *buf = AP_alloc(remoteID,444,sizeof(oldinterfaces_comm));
          oldinterfaces_comm *castbuf = (oldinterfaces_comm *) buf; 
          castbuf->oldnum = myrank;
          castbuf->oldpv  = (*recup)[i].second;
          AP_send(buf);
          sendcounts[remoteID]++;       
        } 
      }
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        oldinterfaces_comm * castbuf = (oldinterfaces_comm *) msg;
        pVertex pv    = castbuf->oldpv;
        int recvID    = castbuf->oldnum;
        int tmp;
        int isChange = EN_getDataInt((pEntity) pv , tagChange, &tmp);
        if(isChange) {
          void *temp_ptr; 
          int isInte = EN_getDataPtr((pEntity) pv , tagNew, &temp_ptr);
          assert(isInte);
          VectorOfCommonpVertex_type *recup = 
            (VectorOfCommonpVertex_type *) temp_ptr;
          unsigned int j=0;
          for( j=0 ; j<(*recup).size() ; j++) {
            int remoteID  = (*recup)[j].first;
            if(remoteID == recvID) break;
          }
          assert(j!= (*recup).size());  
      
        } else {
          void *temp_ptr; 
          int isInte = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
          if(isInte) {
            unsigned int j=0;
            const VectorOfCommonpVertex_type *recup = 
              (const VectorOfCommonpVertex_type *) temp_ptr;
            for( j=0 ; j<(*recup).size() ; j++) {
              int remoteID  = (*recup)[j].first;
              if(remoteID == recvID) break;
            }
            if(j==(*recup).size()) {
              void *temp_ptr2; 
              int is = EN_getDataPtr((pEntity) pv , tagNew, &temp_ptr2);
              assert(is);
              VectorOfCommonpVertex_type *recup2 = 
                (VectorOfCommonpVertex_type *) temp_ptr2;
              if((*recup2).size()==1) break;
              for( j=0 ; j<(*recup2).size() ; j++) {
                int remoteID  = (*recup2)[j].first;
                if(remoteID == recvID) break;
              }
              if(j==(*recup2).size()) {
                for( j=0 ; j<(*recup2).size() ; j++) {
                  printf("recup2 %d\n",(*recup2)[j].first);
                }       
                printf("je suis %d et tu es %d\n",myrank,recvID);
              }  
              assert(j<(*recup2).size());
              int minproc = nproc + 10;
              for( j=0 ; j<(*recup2).size() ; j++) {
                int remoteID  = (*recup2)[j].first;
                minproc = std::min(minproc,remoteID);
              }
              assert(minproc==recvID);
            }
          }
        }
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  }

  // -------------------------------------------------------------------
  void checkNew(pMesh mesh, pMeshDataId tagData,pMeshDataId tagChange ) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
  
    while ((pv = VIter_next(vit))) {
      int tmp;
      int isMaster = EN_getDataInt((pEntity) pv , tagMaster, &tmp);
      if(!isMaster) continue;
      int isChange = EN_getDataInt((pEntity) pv , tagChange, &tmp);
      if(!(isChange && (tmp==1))) continue;
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      assert(isInterface);
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      int sizeinterface = (*recup).size();
      if(sizeinterface<=1) {
        printf("--------------- my %d\n",myrank);
        unsigned int kk;
        for(kk=0 ; kk<(*recup).size() ; kk++) {
          printf("recup %d \n",(*recup)[kk].first);
        }
        void *temp_ptr2; 
        int is = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr2);
        assert(is);
        VectorOfCommonpVertex_type *recup2 = 
          (VectorOfCommonpVertex_type *) temp_ptr2;
        assert((*recup2).size());
        for(kk=0 ; kk<(*recup2).size() ; kk++){
          printf("recup2 %d \n",(*recup2)[kk].first); 
        }
      }
      assert(sizeinterface>1);
      int i;
      for(i=0 ; i<sizeinterface ; i++){
        if(myrank==(*recup)[i].first) break;
      }
      assert(i<sizeinterface);
    }
    VIter_delete(vit); 
  
    vit = M_vertexIter(mesh);  
    while ((pv = VIter_next(vit))) {
      int tmp;
      int isChange = EN_getDataInt((pEntity) pv , tagChange, &tmp);
      if(!(isChange && (tmp==1))) continue;
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      assert(isInterface);
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      int sizeinterface = (*recup).size();
      if(sizeinterface<=1) {
        printf("--------------- my %d\n",myrank);
        unsigned int kk;
        for(kk=0 ; kk<(*recup).size() ; kk++){
          printf("recup %d \n",(*recup)[kk].first);
        }
        void *temp_ptr2; 
        int is = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr2);
        assert(is);
        VectorOfCommonpVertex_type *recup2 =
          (VectorOfCommonpVertex_type *) temp_ptr2;
        assert((*recup2).size());
        for(kk=0 ; kk<(*recup2).size() ; kk++) {
          printf("recup2 %d %p\n",(*recup2)[kk].first,(*recup2)[kk].second);   
        } 
      } 
      assert(sizeinterface>1);
      int i;
      for(i=0 ; i<sizeinterface ; i++){
        if(myrank==(*recup)[i].first) break;
      }
      assert(i<sizeinterface);
    }
    VIter_delete(vit);   
  }

#endif

  // -------------------------------------------------------------------
  void UpdateInterfaces4(pMesh mesh) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isMaster = EN_getDataInt((pEntity) pv , tagMaster,&tmp);
      if(!isMaster) continue;
    
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      assert(tmp==1 && isChanged);

      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      assert(isInterface);
      const VectorOfCommonpVertex_type *recup =
        (const VectorOfCommonpVertex_type *) temp_ptr;
      assert((*recup).size()>1);

      for(unsigned int i=0 ; i<(*recup).size() ; i++) {
        int remoteID = (*recup)[i].first;
        if(remoteID==myrank) continue;
        pVertex remoteP = (*recup)[i].second;
        assert(remoteP);
        for(unsigned int j=0 ; j<(*recup).size() ; j++) {
          int intID = (*recup)[j].first;
          if(intID==remoteID || intID==myrank) continue;
          pVertex intP = (*recup)[j].second;
          assert(intP);
          void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
          points_comm_new *castbuf = (points_comm_new *) buf; 
          castbuf->pdest  = remoteP;
          castbuf->psend  = intP;
          castbuf->sendID = intID;
          AP_send(buf);
          sendcounts[remoteID]++;      
        }
      }
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        points_comm_new * castbuf = (points_comm_new *) msg;
        pVertex pmine =  castbuf->pdest;
        pVertex precv  = castbuf->psend;
        int     recvID = castbuf->sendID;
        void *temp_ptr; 
        int isInterface = EN_getDataPtr((pEntity) pmine , tagData, &temp_ptr);
        int tmp; 
        int isMaster = EN_getDataInt((pEntity) pmine , tagMaster,&tmp);
        assert(isInterface);
        assert(!isMaster);
        VectorOfCommonpVertex_type *recup = 
          (VectorOfCommonpVertex_type *) temp_ptr;
        unsigned int j=0;
        for(j=0 ; j<(*recup).size() ; j++) {
          int remoteID  = (*recup)[j].first;
          if(remoteID == recvID) {
            (*recup)[j].second = precv;
            break;
          }
        }           
        assert(j<(*recup).size());    
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  
    vit = M_vertexIter(mesh);  
    while ((pv = VIter_next(vit))) {
      int tmp;     
      EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(tmp!=1) continue;

      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      assert(isInterface);
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      assert((*recup).size()>1);

      VectorOfCommonpVertex_type newInterface;
      for(unsigned int i=0 ; i<(*recup).size() ; i++) {
        int remoteID = (*recup)[i].first;
        if(remoteID==myrank) continue;
        pVertex remoteP = (*recup)[i].second;
        assert(remoteP);
        newInterface.push_back(CommonpVertex_type(remoteID,remoteP));
      }
      EN_attachDataPtr((pEntity) pv ,tagData, 
                       new VectorOfCommonpVertex_type(newInterface));                

    }
    VIter_delete(vit);   


  }  
  
  // -------------------------------------------------------------------
  void UpdateInterfaces3(pMesh mesh) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(!isChanged) continue;
      if(tmp==2) continue;
      if(tmp==3) continue;
      assert(tmp==1);
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      int isMaster = EN_getDataInt((pEntity) pv , tagMaster,&tmp);
      if(isMaster) assert(isInterface);
      assert(isInterface);
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      if(isMaster) assert((*recup).size()>1);
      if(!isMaster) {
        assert((*recup).size()>1);
        int remoteID = nproc + 10;
        pVertex psend;
        for(unsigned int i=0 ; i<(*recup).size() ; i++) {
          if(remoteID > (*recup)[i].first) {
            remoteID = (*recup)[i].first;
            psend    = (*recup)[i].second;
          }
        }
        assert(remoteID < nproc);
        assert(remoteID!=myrank);
        assert(remoteID<myrank);
        for(unsigned int i=0 ; i<(*recup).size() ; i++) {
          assert(remoteID <=(*recup)[i].first);
        }
        void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
        points_comm_new *castbuf = (points_comm_new *) buf; 
        castbuf->pdest  = psend;
        castbuf->psend  = pv;
        castbuf->sendID = myrank;
        AP_send(buf);
        sendcounts[remoteID]++;      
      }
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        points_comm_new * castbuf = (points_comm_new *) msg;
        pVertex pmine =  castbuf->pdest;
        pVertex precv  = castbuf->psend;
        int     recvID = castbuf->sendID;
        void *temp_ptr; 
        int isInterface = EN_getDataPtr((pEntity) pmine , tagData, &temp_ptr);
        int tmp; 
        int isMaster = EN_getDataInt((pEntity) pmine , tagMaster,&tmp);
        assert(isInterface);
        assert(isInterface && isMaster);
        VectorOfCommonpVertex_type *recup = 
          (VectorOfCommonpVertex_type *) temp_ptr;
        unsigned int j=0;
        for(j=0 ; j<(*recup).size() ; j++) {
          int remoteID  = (*recup)[j].first;
          if(remoteID == recvID) {
            (*recup)[j].second = precv;
            break;
          }
        }           
        assert(j<(*recup).size());    
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  }  
  
  // -------------------------------------------------------------------
  void UpdateInterfaces2(pMesh mesh, MDB_DataExchanger &de) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
 
    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isMaster = EN_getDataInt((pEntity) pv , tagMaster,&tmp);
      if(!isMaster) continue;
#ifdef DEBUG
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      assert(isChanged);
      assert(tmp==1);
#endif
      void *temp_ptr2; 
      int is = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr2);
      assert(is);
      const VectorOfCommonpVertex_type *recup2 = 
        (const VectorOfCommonpVertex_type *) temp_ptr2;
      int sizeinterface = (*recup2).size();
      assert(sizeinterface>1);
      int  *listID=new int[sizeinterface];
      for(int i=0 ; i<sizeinterface ; i++){
        listID[i] = (*recup2)[i].first;
      }
#ifdef DEBUG
      int k;
      for(k=0 ; k<sizeinterface ; k++){
        if(myrank==(*recup2)[k].first) break;
      }
      assert(k!=sizeinterface);
#endif
      void *temp_ptr; 
      int oldPoint = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr);
      VectorOfCommonpVertex_type *recup = 
        (VectorOfCommonpVertex_type *) temp_ptr;
      assert(oldPoint);
      assert((*recup).size());
#ifdef DEBUG
      unsigned int kk;
      for(kk=0 ; kk<(*recup).size() ; kk++){
        if(myrank==(*recup)[kk].first) break;
      }
      if(kk!=(*recup).size()) {
        if((*recup).size()!=1) {
          printf("----------------- myrank %d %d\n",myrank,kk);
          for(k=0 ; k<sizeinterface ; k++){
            printf("recup2 %d \n",(*recup2)[k].first);
          }
          for(kk=0 ; kk<(*recup).size() ; kk++){
            printf("recup %d \n",(*recup)[kk].first);
          }
        }
        assert((*recup).size()==1);
      }
#endif

      int *tabp=new int[nproc];
      int *tabproc=new int[nproc];
      for(int i=0 ; i<nproc ; i++) tabproc[i] = 0;
      for(unsigned int i=0 ; i<(*recup).size() ; i++) {
        tabproc[(*recup)[i].first] += 1;
        tabp[(*recup)[i].first]     = i;
      }
      for(int i=0 ; i<sizeinterface ; i++) {
        int remoteID = (*recup2)[i].first;
        if(remoteID==myrank) continue;
        if(tabproc[remoteID]) { /*already exists --  send my pointer, my ID and its pointer*/
          assert(tabproc[remoteID]==1);
          assert(remoteID!=myrank);
          int sizebuf;
          void *msg = de.sendData ((pEntity) pv, remoteID, sizebuf );
  
          void *buf = AP_alloc(remoteID,de.tag(),sizeof(int) + 
                               sizeof(points_comm_new) + sizebuf);
          char *castbuf = (char *) buf; 
          int  cas = -1;
          memcpy(&castbuf[0],&cas,sizeof(int));
          points_comm_new pcn;
          pcn.pdest  = ((*recup)[tabp[remoteID]]).second;
          assert(((*recup)[tabp[remoteID]]).first==remoteID);
          pcn.psend  = pv;
          pcn.sendID = myrank;
          pcn.destID = remoteID;
          memcpy(&castbuf[sizeof(int)],&pcn,sizeof(points_comm_new));
          memcpy(&castbuf[sizeof(int)+sizeof(points_comm_new)],msg,sizebuf);        
          free(msg);
          AP_send(buf);
          sendcounts[remoteID]++;       
        } else { /*new for proc remoteID  -- send coor point + Remote info + my pointer*/
          assert(remoteID!=myrank);
          int sizebuf;
          void *msg = de.sendData ((pEntity) pv, remoteID, sizebuf );
          void *buf = AP_alloc(remoteID,de.tag(),2*sizeof(int)+ 
                               sizeof(coor_comm)+ (sizeinterface) * sizeof(int) + sizebuf );
          char *castbuf = (char *) buf; 
          int  cas = 1;
          memcpy(&castbuf[0],&cas,sizeof(int));
          coor_comm  coorcom;           
          coorcom.X     = P_x(pv);
          coorcom.Y     = P_y(pv);
          coorcom.Z     = P_z(pv);
          pGEntity  pg  = EN_whatIn(pv);
          coorcom.tag   = GEN_tag(pg);
          coorcom.dim   = GEN_type(pg);
          coorcom.nproc = myrank;
          coorcom.psend = pv;
          memcpy(&castbuf[sizeof(int)],&coorcom,sizeof(coor_comm));
          memcpy(&castbuf[sizeof(int) + sizeof(coor_comm)],&sizeinterface,sizeof(int));
          memcpy(&castbuf[sizeof(int) + sizeof(coor_comm) + sizeof(int)],&listID[0],sizeinterface * sizeof(int));
          memcpy(&castbuf[2*sizeof(int) + sizeof(coor_comm)+ (sizeinterface) * sizeof(int)],msg,sizebuf);
          free(msg);
          AP_send(buf);
          sendcounts[remoteID]++;                  
        }
      }       
      delete []tabp;
      delete []tabproc;
      delete []listID;
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        char * castbuf = (char *) msg;
        int cas = *((int *) &castbuf[0]);
        if(cas == 1) {/*point doesn't exist*/
          coor_comm * coorcom = (coor_comm * ) &castbuf[sizeof(int)];
          pGEntity pent;
          int dim  = coorcom->dim;
          int tag  = coorcom->tag;
          if(dim == 0) {
            pent = (pGEntity) GM_vertexByTag(mesh->model,tag);
          } else if(dim==1) {
            pent = (pGEntity) GM_edgeByTag(mesh->model,tag);     
          } else if(dim==2) {
            pent = (pGEntity) GM_faceByTag(mesh->model,tag);
          } else if(dim==3) {
            pent = (pGEntity) GM_regionByTag(mesh->model,tag);
          } else {
            printf("pbs**** %d\n",dim);
          }
          //           double  X = coorcom->X;
          //           double  Y = coorcom->Y;
          //           double  Z = coorcom->Z;
          //           pVertex pnew = M_createVP(mesh,X,Y,Z,-1,pent);
          
          double  XYZ[3] = {coorcom->X,coorcom->Y,coorcom->Z};
          pVertex pnew = M_createV2(mesh,XYZ,-1,pent);

          assert(pnew);
          int procdep  = coorcom->nproc;
          pVertex pdep = coorcom->psend;
          unsigned int sizeinterfaces = *((int *) &castbuf[sizeof(int) + 
                                                           sizeof(coor_comm)]);
          assert(sizeinterfaces > 1);
          VectorOfCommonpVertex_type RemoteData;
          int isOK=0;
          for(unsigned int i=0 ; i<sizeinterfaces ; i++) {
            int num = *((int *) &castbuf[2*sizeof(int) + sizeof(coor_comm) 
                                         + i*sizeof(int)]);
            if(num==procdep) {
              isOK = 1;
              RemoteData.push_back(CommonpVertex_type(procdep,pdep));
            } else {
              RemoteData.push_back(CommonpVertex_type(num,NULL));     
            }
          } 
          assert(isOK);
          EN_attachDataInt((pEntity) pnew ,tagChange, 1);
          assert(RemoteData.size()==sizeinterfaces);
          assert(RemoteData.size()>1);

          EN_attachDataPtr((pEntity) pnew ,tagData, 
                           new VectorOfCommonpVertex_type(RemoteData));               
          de.receiveData (pnew,from, &castbuf[2*sizeof(int)
                                              + sizeof(coor_comm)+ (sizeinterfaces) * sizeof(int)]);
           
        } else { /*point already exists*/
          assert(cas==-1);
          points_comm_new * pcn = (points_comm_new * ) &castbuf[sizeof(int)];
          pVertex pv    = pcn->pdest;
          pVertex precv = pcn->psend;
          int recvID    = pcn->sendID;
          int destID    = pcn->destID;
          assert(pv);
          assert(destID==myrank);
          assert(recvID!=myrank);
          void *temp_ptr; 
          int isInte = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr);
          assert(isInte);
#ifdef DEBUG
          int tmp; 
          int isMaster = EN_getDataInt((pEntity) pv , tagMaster,&tmp);
          assert(isInte && !isMaster);
          int isC = EN_getDataInt((pEntity) pv , tagChange,&tmp);
          assert(isC && (tmp==1));
#endif
          VectorOfCommonpVertex_type *recup = 
            (VectorOfCommonpVertex_type *) temp_ptr;
          assert((*recup).size()>1);
          VectorOfCommonpVertex_type newRemote;
          for(unsigned int j=0 ; j<(*recup).size() ; j++) {
            int remoteID  = (*recup)[j].first;
            if(remoteID == recvID) {
              newRemote.push_back(CommonpVertex_type(remoteID,precv));
            } else {
              newRemote.push_back(CommonpVertex_type(remoteID,NULL));
            }
          }
          assert(newRemote.size()>1);
          void *temp_ptr2; 
          int is = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr2);
          if(is) EN_deleteData((pEntity) pv , tagData);

          EN_attachDataPtr((pEntity) pv ,tagData, 
                           new VectorOfCommonpVertex_type(newRemote));                
          de.receiveData (pv,from, &castbuf[sizeof(int)+sizeof(points_comm_new)]);
          
        }     
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  
  }

  
  // -------------------------------------------------------------------
  void UpdateInterfaces1(pMesh mesh, MDB_DataExchanger &de)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
 
    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(!isChanged) continue;
      void *temp_ptr2; 
      int is = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr2);
      assert(is);
      VectorOfCommonpVertex_type *recup2 = 
        (VectorOfCommonpVertex_type *) temp_ptr2;
      if((*recup2).size()==1) {
   
        EN_attachDataInt((pEntity) pv , tagChange,2);
        continue;
      }
      assert((*recup2).size()>1);
      /*est ce que je suis a la fin?*/
      unsigned int i=0;
      for(i=0 ; i< (*recup2).size(); i++){
        if(myrank==(*recup2)[i].first) break;
      }
      if(i==(*recup2).size()) EN_attachDataInt((pEntity) pv , tagChange,3);
     
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
    
      int remoteID;
    
      if(isInterface) {
        const VectorOfCommonpVertex_type *recup = 
          (const VectorOfCommonpVertex_type *) temp_ptr;
        int minproc = nproc + 10;
        assert((*recup).size());
        for(unsigned int i=0 ; i<(*recup).size() ; i++){
          minproc =         std::min(minproc,(*recup)[i].first);
        }
        minproc = std::min(minproc,myrank);
       
        remoteID = nproc + 10 ;
        for(unsigned int i=0 ; i<(*recup2).size() ; i++) {
          remoteID =        std::min(remoteID,(*recup2)[i].first);
        }
        assert(remoteID<nproc && remoteID > -1);
      
        if(remoteID == myrank) {
          EN_attachDataInt((pEntity) pv , tagMaster,1);
          /*swap RemoteData and Data*/
          VectorOfCommonpVertex_type temp ;
          for(unsigned int i=0 ; i<(*recup).size() ; i++) {
            int   num = (*recup)[i].first;
            pVertex p = (*recup)[i].second;
            temp.push_back(CommonpVertex_type(num,p));
          }
          assert((*recup2).size()>1);
          VectorOfCommonpVertex_type temp2 ;
          for(unsigned int i=0 ; i<(*recup2).size() ; i++) {
            int   num = (*recup2)[i].first;
            pVertex p = (*recup2)[i].second;
            temp2.push_back(CommonpVertex_type(num,p));
          }
          delete recup2;
          delete recup;
          EN_attachDataPtr((pEntity) pv , tagData, 
                           new VectorOfCommonpVertex_type(temp2));    
          EN_attachDataPtr((pEntity) pv , tagRemote, 
                           new VectorOfCommonpVertex_type(temp));           
          continue;
        }
      
        if(minproc != myrank) continue; 
        assert(minproc==myrank);
        int recupsize  = (*recup).size();
        assert(recupsize);
        int recupsize2 = (*recup2).size();
        int i=0;
        for(i=0 ; i< recupsize; i++){
          if(remoteID==(*recup)[i].first) break;
        }
        if(i!=recupsize) continue;       /*proc remoteID already knows info*/
        assert(i==recupsize);
        /*send info to remoteID*/
        int sizebuf;
        void *msg = de.sendData ((pEntity) pv, remoteID, sizebuf );
        void *buf = AP_alloc(remoteID,de.tag(),2*sizeof(int)+ recupsize 
                             * sizeof(oldinterfaces_comm)
                             + recupsize2 * sizeof(newinterfaces_comm)
                             +  sizeof(coor_comm) + sizebuf);
        char *castbuf = (char *) buf; 
        memcpy(&castbuf[0],&(recupsize),sizeof(int));          
        for(i=0 ; i<recupsize ; i++) {  
          oldinterfaces_comm oldint;
          oldint.oldnum = (*recup)[i].first;
          oldint.oldpv  = (*recup)[i].second;          
          memcpy(&castbuf[i*sizeof(oldinterfaces_comm) 
                          + sizeof(int)],&oldint,sizeof(oldinterfaces_comm));
        }
        memcpy(&castbuf[recupsize*sizeof(oldinterfaces_comm) 
                        + sizeof(int)],&(recupsize2),sizeof(int));           
        for(i=0 ; i<recupsize2 ; i++) {
          newinterfaces_comm newint;
          newint.newproc = (*recup2)[i].first;
          memcpy(&castbuf[i * sizeof(newinterfaces_comm) 
                          + recupsize*sizeof(oldinterfaces_comm) + 2*sizeof(int)],
                 &newint,sizeof(newinterfaces_comm));
        }
        coor_comm  coorcom;       
        coorcom.X     = P_x(pv);
        coorcom.Y     = P_y(pv);
        coorcom.Z     = P_z(pv);
        pGEntity  pg  = EN_whatIn(pv);
        coorcom.tag   = GEN_tag(pg);
        coorcom.dim   = GEN_type(pg);
        coorcom.nproc = myrank;
        coorcom.psend = pv;
        memcpy(&castbuf[recupsize * sizeof(oldinterfaces_comm)
                        + recupsize2 * sizeof(newinterfaces_comm) 
                        + 2*sizeof(int)],&coorcom,sizeof(coor_comm));                    
        memcpy(&castbuf[recupsize * sizeof(oldinterfaces_comm)
                        + recupsize2 * sizeof(newinterfaces_comm) + 2*sizeof(int)
                        + sizeof(coor_comm)],msg,sizebuf);                     
        free(msg);
        AP_send(buf);
        sendcounts[remoteID]++; 
      } else {    
        remoteID = nproc + 10;
        assert((*recup2).size() != 1);
        for(unsigned int i=0 ; i<(*recup2).size() ; i++){ 
          remoteID = std::min(remoteID,(*recup2)[i].first);
        }
        assert(remoteID<nproc);
        if(remoteID == myrank) {
          int tmp;
          int is = EN_getDataInt((pEntity) pv , tagMaster,&tmp);
          assert(!is);
          EN_attachDataInt((pEntity) pv , tagMaster, 1);
          /*swap RemoteData and Data*/
          assert((*recup2).size()>1);
          VectorOfCommonpVertex_type temp2 ;
          for(unsigned int i=0 ; i<(*recup2).size() ; i++) {
            int   num = (*recup2)[i].first;
            pVertex p = (*recup2)[i].second;
            temp2.push_back(CommonpVertex_type(num,p));
          }
          EN_attachDataPtr((pEntity) pv , tagData, 
                           new VectorOfCommonpVertex_type(temp2));    
          VectorOfCommonpVertex_type temp;
          temp.push_back(CommonpVertex_type(myrank,pv));
          EN_attachDataPtr((pEntity) pv , tagRemote, 
                           new VectorOfCommonpVertex_type(temp));     
           
  
          continue;
        }   
        /*send info to remoteID*/
        int recupsize = 0;
        int recupsize2 = (*recup2).size();
        int sizebuf;
        void *msg = de.sendData ((pEntity) pv, remoteID, sizebuf );
        void *buf = AP_alloc(remoteID,de.tag(),2*sizeof(int)
                             + recupsize2 * sizeof(newinterfaces_comm)
                             +  sizeof(coor_comm) +sizebuf);
        char *castbuf = (char *) buf; 
        memcpy(&castbuf[0],&(recupsize),sizeof(int));          
        memcpy(&castbuf[sizeof(int)],&(recupsize2),sizeof(int));           
        for(int i=0 ; i<recupsize2 ; i++)            
          memcpy(&castbuf[i * sizeof(newinterfaces_comm)  + 2*sizeof(int) ],
                 &(*recup2)[i].first,sizeof(newinterfaces_comm));
        coor_comm  coorcom;       
        coorcom.X     = P_x(pv);
        coorcom.Y     = P_y(pv);
        coorcom.Z     = P_z(pv);
        pGEntity  pg  = EN_whatIn(pv);
        coorcom.tag   = GEN_tag(pg);
        coorcom.dim   = GEN_type(pg);
        coorcom.nproc = myrank;
        coorcom.psend = pv;
        memcpy(&castbuf[recupsize * sizeof(oldinterfaces_comm)
                        + recupsize2 * sizeof(newinterfaces_comm) 
                        + 2*sizeof(int) ],&coorcom,sizeof(coor_comm));                     
        memcpy(&castbuf[2*sizeof(int)+  recupsize2 * sizeof(newinterfaces_comm)
                        + sizeof(coor_comm)],msg,sizebuf);                     
        free(msg);
        AP_send(buf);
        sendcounts[remoteID]++; 
      }
    }

    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        char * castbuf = (char *) msg;
        VectorOfCommonpVertex_type oldRemoteData;
        int sizeold = *((int *) &castbuf[0]);
        for(int i=0 ; i<sizeold ; i++) {
          oldinterfaces_comm * oldint = 
            (oldinterfaces_comm *) &castbuf[sizeof(int)+ i*sizeof(oldinterfaces_comm)];
          int    num = oldint->oldnum;
          pVertex pv = oldint->oldpv;
          oldRemoteData.push_back(CommonpVertex_type(num,pv));
        }      
        VectorOfCommonpVertex_type newRemoteData;
        int sizenew = *((int *) &castbuf[sizeof(int) + sizeold 
                                         * sizeof(oldinterfaces_comm)]);
        assert(sizenew>1);
        for(int i=0 ; i<sizenew ; i++) {
          newinterfaces_comm * newint = 
            (newinterfaces_comm *) &castbuf[2*sizeof(int)+ sizeold 
                                            * sizeof(oldinterfaces_comm) + i * sizeof(newinterfaces_comm)];
          int nn = newint->newproc; 
          newRemoteData.push_back(CommonpVertex_type(nn,NULL));
      
        }
        coor_comm * coorcom = (coor_comm * ) &castbuf[2*sizeof(int) 
                                                      +sizeold*sizeof(oldinterfaces_comm) + 
                                                      sizenew*sizeof(newinterfaces_comm)];
        pGEntity pent;
        int dim  = coorcom->dim;
        int tag  = coorcom->tag;
        if(dim == 0) {
          pent = (pGEntity) GM_vertexByTag(mesh->model,tag);
        } else if(dim==1) {
          pent = (pGEntity) GM_edgeByTag(mesh->model,tag);     
        } else if(dim==2) {
          pent = (pGEntity) GM_faceByTag(mesh->model,tag);
        } else if(dim==3) {
          pent = (pGEntity) GM_regionByTag(mesh->model,tag);
        } else {
          printf("pbs**** %d\n",dim);
          pent = (pGEntity) GM_vertexByTag(mesh->model,tag);
        }
        
        //         double  X = coorcom->X;
        //         double  Y = coorcom->Y;
        //         double  Z = coorcom->Z;
        //         pVertex pnew = M_createVP(mesh,X,Y,Z,-1,pent);
        
        double  XYZ[3] = {coorcom->X,coorcom->Y,coorcom->Z};
        pVertex pnew = M_createV2(mesh,XYZ,-1,pent);
        assert(pnew);
        int procdep = coorcom->nproc;
        //  if(sizeold!=0) assert(procdep==-1);
        //   if(sizeold == 0) { 
        //  assert(!oldRemoteData.size());     
        pVertex pinit = coorcom->psend;
        oldRemoteData.push_back(CommonpVertex_type(procdep,pinit));
        //  }

        EN_attachDataInt((pEntity) pnew , tagChange,1);
        EN_attachDataInt((pEntity) pnew , tagMaster,1);
        assert(newRemoteData.size()>1);
       
        EN_attachDataPtr((pEntity) pnew , tagData, 
                         new VectorOfCommonpVertex_type(newRemoteData));
        assert(oldRemoteData.size());     
        EN_attachDataPtr((pEntity) pnew , tagRemote, 
                         new VectorOfCommonpVertex_type(oldRemoteData));    
            
        de.receiveData (pnew,from, &castbuf[2*sizeof(int)+ sizeold * 
                                            sizeof(oldinterfaces_comm)
                                            + sizenew * sizeof(newinterfaces_comm)
                                            + sizeof(coor_comm)]);

        AP_free(msg);
      }    
    }  
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;

  }

  // -------------------------------------------------------------------
  void UpdateInterfaces(pMesh mesh, MDB_DataExchanger &de) {

#ifdef DEBUG
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
#endif

    /*Phase 1 : send info to the smallest rank*/
    UpdateInterfaces1(mesh,de);
#ifdef DEBUG
    checkRemotePointerChange(mesh,tagData,tagRemote,tagMaster);
    puts("check third ok");
#endif

    /*Phase 2 : send info to all proc*/
    UpdateInterfaces2(mesh,de);
#ifdef DEBUG  
    checkNew(mesh,tagData,tagChange);
    puts("check forth ok");
#endif  

    /*Phase 3 : update info of smallest rank*/
    UpdateInterfaces3(mesh);
  
    /*Phase 4 : update info of all proc*/
    UpdateInterfaces4(mesh);
  }


  // -------------------------------------------------------------------
  void UpdatePointerInterface(pMesh mesh) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(!isChanged) continue;
      if((tmp==2) || (tmp==10) || (tmp==11)) continue; //if 1 : interface new -- if 3 : must be deleted
      assert(tmp==1 || tmp==3);
      int isMaster = EN_getDataInt((pEntity) pv , tagMaster,&tmp);
      if(!isMaster) continue;

      void *temp_ptr2; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr2);
      assert(isInterface);
      const VectorOfCommonpVertex_type *recup2 = 
        (const VectorOfCommonpVertex_type *) temp_ptr2;
      assert((*recup2).size());
    
      void *temp_ptr; 
      int isInternal = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr);
      assert(isInternal);
      VectorOfCommonpVertex_type *recup = (VectorOfCommonpVertex_type *) temp_ptr;
      int sizeold = (*recup).size();
      if(sizeold==1){//old internal point
        int remoteID = (*recup)[0].first;
        if(remoteID==myrank) continue;
        pVertex remoteP  = (*recup)[0].second;
        assert(remoteP);
        for(unsigned int i=0 ; i <(*recup2).size() ; i++) {
          int newID    = (*recup2)[i].first;
          pVertex newP = (*recup2)[i].second;
          assert(newP);
          void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
          points_comm_new *castbuf = (points_comm_new *) buf; 
          castbuf->pdest  = remoteP;
          castbuf->destID = remoteID;
          castbuf->psend  = newP;
          castbuf->sendID = newID;
          AP_send(buf);
          sendcounts[remoteID]++;     
        }
        void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
        points_comm_new *castbuf = (points_comm_new *) buf; 
        castbuf->pdest  = remoteP;
        castbuf->destID = remoteID;
        castbuf->psend  = pv;
        castbuf->sendID = myrank;
        AP_send(buf);
        sendcounts[remoteID]++;       
      
      } else { //old interface point
        int *tabproc=new int[nproc];
        for(int i=0 ; i<nproc ; i++) tabproc[i] = 0;
        for(unsigned int i=0 ; i <(*recup2).size() ; i++) {
          tabproc[(*recup2)[i].first] = 1;
        }
        for(unsigned int i=0 ; i <(*recup).size() ; i++) {
          int remoteID = (*recup)[i].first;
          if(tabproc[remoteID]) continue;
          pVertex remoteP = (*recup)[i].second;
          assert(remoteP);
          for(unsigned int j=0 ; j <(*recup2).size() ; j++) {
            int newID    = (*recup2)[j].first;
            pVertex newP = (*recup2)[j].second;
            assert(newID!=remoteID); assert(newID!=myrank);
            assert(newP);
            void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
            points_comm_new *castbuf = (points_comm_new *) buf; 
            castbuf->pdest  = remoteP;
            castbuf->destID = remoteID;
            castbuf->psend  = newP;
            castbuf->sendID = newID;
            AP_send(buf);
            sendcounts[remoteID]++;     
          }
          void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
          points_comm_new *castbuf = (points_comm_new *) buf; 
          castbuf->pdest  = remoteP;
          castbuf->destID = remoteID;
          castbuf->psend  = pv;
          castbuf->sendID = myrank;
          AP_send(buf);
          sendcounts[remoteID]++;         
        }
        delete []tabproc;
      }
     
    }
    VIter_delete(vit);   
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        points_comm_new * pcn = (points_comm_new * ) msg;
        pVertex precv = pcn->pdest;
        assert(myrank == pcn->destID);
        int procdep = pcn->sendID;
        pVertex pv  = pcn->psend;
        void *temp_ptr; 
        int isInternal = EN_getDataPtr((pEntity) precv , tagRemote, &temp_ptr);
        assert(isInternal);
        VectorOfCommonpVertex_type *recup = (VectorOfCommonpVertex_type *) temp_ptr;
#ifdef DEBUG
        int tmp;
        int isMaster = EN_getDataInt((pEntity) precv , tagMaster, &tmp);
        assert(!isMaster);
#endif
        unsigned int i;
        for(i=0 ; i<(*recup).size() ; i++) {
          int remoteID = (*recup)[i].first;
          if(remoteID == procdep) {
            (*recup)[i].second = pv;
            break;
          }      
        }
        assert(i<(*recup).size());    
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  
  }
  // -------------------------------------------------------------------
  void UpdateoldInterface(pMesh mesh) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    //   pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");

    int *  sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(!isChanged) continue;
      if(!(tmp==2 || tmp==10 || tmp==11)) continue; //if 1 : interface new -- if 3 : must be deleted
      void *temp_ptr2; 
      int is = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr2);
      assert(is);
      VectorOfCommonpVertex_type *recup2 = (VectorOfCommonpVertex_type *) temp_ptr2;
      assert((*recup2).size()==1);
    
      void *temp_ptr; 
      int isoldInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      if(!isoldInterface) continue;
      const VectorOfCommonpVertex_type *recup = 
        (const VectorOfCommonpVertex_type *) temp_ptr;
      assert((*recup).size());
    
      int minproc = nproc + 10;
      for(unsigned int i=0 ; i<(*recup).size() ; i++){
        minproc = std::min(minproc,(*recup)[i].first);
      }
      minproc = std::min(minproc,myrank);
      if(minproc!=myrank) continue;
      int newID = (*recup2)[0].first;
      pVertex newP = (*recup2)[0].second;
      if(!newP) {
        assert(newID==myrank);
        newP = pv;
      }
    
      for(unsigned i=0 ; i<(*recup).size() ; i++) {
        int remoteID = (*recup)[i].first;
        assert(remoteID!=myrank);
        pVertex remoteP = (*recup)[i].second;
        void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
        points_comm_new *castbuf = (points_comm_new *) buf; 
        castbuf->pdest  = remoteP;
        castbuf->destID = remoteID;
        castbuf->psend  = newP;
        castbuf->sendID = newID;
        AP_send(buf);
        sendcounts[remoteID]++;         
      }
    }
    VIter_delete(vit);   
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        points_comm_new * pcn = (points_comm_new * ) msg;
        pVertex precv = pcn->pdest;
        assert(myrank == pcn->destID);
        int procdep = pcn->sendID;
        pVertex pv  = pcn->psend;
        void *temp_ptr; 
        int isInternal = EN_getDataPtr((pEntity) precv , tagRemote, &temp_ptr);
        assert(isInternal);
        VectorOfCommonpVertex_type *recup = (VectorOfCommonpVertex_type *) temp_ptr;
#ifdef DEBUG
        int tmp;
        int isChange = EN_getDataInt((pEntity) precv , tagChange, &tmp);
        assert(isChange);
        //assert(tmp==2);
        assert((*recup).size()==1);
#endif
        int remoteID = (*recup)[0].first;
        assert(remoteID == procdep);
        (*recup)[0].second = pv;
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;  


  }

  // -------------------------------------------------------------------
  void UpdatePointerInternal(pMesh mesh) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    //   pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(!isChanged) continue;
      if(!(tmp==10 || tmp==11)) continue; //if 1 : interface new -- if 3 : must be deleted

      void *temp_ptr2; 
      int is = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr2);
      assert(is);
      VectorOfCommonpVertex_type *recup2 = (VectorOfCommonpVertex_type *) temp_ptr2;
      assert((*recup2).size()==1);
    
      if(tmp==10) {
        int remoteID = (*recup2)[0].first;
        if(remoteID==myrank) continue;
        pVertex remoteP = (*recup2)[0].second;
        assert(remoteP);
        void *buf = AP_alloc(remoteID,444,sizeof(points_comm_new));
        points_comm_new *castbuf = (points_comm_new *) buf; 
        castbuf->pdest  = remoteP;
        castbuf->destID = remoteID;
        castbuf->psend  = pv;
        castbuf->sendID = myrank;
        AP_send(buf);
        sendcounts[remoteID]++;
      } else {
        void *temp_ptr; 
        is = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
        assert(is);
        const VectorOfCommonpVertex_type *recup = 
          (const VectorOfCommonpVertex_type *) temp_ptr;
        assert((*recup).size());
        int minproc = nproc + 10;
        int ip ;
        for(unsigned int i=0 ; i<(*recup).size() ; i++) {
          if(minproc > (*recup)[i].first) {
            minproc = (*recup)[i].first;
            ip      = i;
          }
        }   
        minproc = std::min(minproc,myrank);
        if(minproc==myrank) continue;
        void *buf = AP_alloc(minproc,444,sizeof(points_comm_new));
        points_comm_new *castbuf = (points_comm_new *) buf; 
        castbuf->pdest  = (*recup)[ip].second;
        castbuf->destID =  minproc;
        castbuf->psend  = pv;
        castbuf->sendID = myrank;
        AP_send(buf);
        sendcounts[minproc]++;
      }     
    }
    VIter_delete(vit);   
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, 444, AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        points_comm_new * pcn = (points_comm_new * ) msg;
        pVertex precv = pcn->pdest;
        assert(myrank == pcn->destID);
        int procdep = pcn->sendID;
        pVertex pv  = pcn->psend;
        void *temp_ptr; 
        int isInternal = EN_getDataPtr((pEntity) precv , tagRemote, &temp_ptr);
        assert(isInternal);
        VectorOfCommonpVertex_type *recup = (VectorOfCommonpVertex_type *) temp_ptr;
#ifdef DEBUG
        int tmp;
        int isChange = EN_getDataInt((pEntity) precv , tagChange, &tmp);
        assert(isChange);
        assert(tmp==2);
        assert((*recup).size()==1);
#endif
        int remoteID = (*recup)[0].first;
        assert(remoteID == procdep);
        (*recup)[0].second = pv;
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  
    UpdateoldInterface(mesh);
  }

  // -------------------------------------------------------------------
  void SendVertex(pMesh mesh, MDB_DataExchanger &de)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");

    VIter vit = M_vertexIter(mesh);  
    pVertex pv;
    while ((pv = VIter_next(vit))) {
      int tmp; 
      int isChanged = EN_getDataInt((pEntity) pv , tagChange,&tmp);
      if(!isChanged) continue;
      if(tmp!=2) continue; //if 1 : interface new -- if 3 : must be deleted

      void *temp_ptr2; 
      int isInternal = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr2);
      assert(isInternal);
      VectorOfCommonpVertex_type *recup2 = 
        (VectorOfCommonpVertex_type *) temp_ptr2;
      assert((*recup2).size()==1);
        
      void *temp_ptr; 
      int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
      if(isInterface){
        /*send only if I'm minproc*/
        const VectorOfCommonpVertex_type *recup = 
          (const VectorOfCommonpVertex_type *) temp_ptr;
        int minproc = nproc + 10;
        int *tabproc=new int[nproc];
        for(int i=0 ; i<nproc ; i++) tabproc[i] = 0;
        for(unsigned int i=0 ; i<(*recup).size() ; i++) {
          minproc = std::min(minproc,(*recup)[i].first);
          tabproc[(*recup)[i].first] = 1;
        }
        minproc = std::min(minproc,myrank);
      
        int remoteID = (*recup2)[0].first;
        if(remoteID==myrank)  {
          EN_attachDataInt((pEntity) pv , tagChange,11);
          continue;
        }      
      
        if(minproc!=myrank) continue;
      
      
        if(tabproc[remoteID]) continue; //remoteID already knows info
        int sizebuf;
        void *msg = de.sendData ((pEntity) pv, remoteID, sizebuf );
        void *buf = AP_alloc(remoteID,de.tag(),sizeof(coor_comm)+sizebuf);

        char *cast = (char *) buf;  
        coor_comm castbuf; 
        castbuf.X     = P_x(pv);
        castbuf.Y     = P_y(pv);
        castbuf.Z     = P_z(pv);
        pGEntity  pg  = EN_whatIn(pv);
        castbuf.tag   = GEN_tag(pg);
        castbuf.dim   = GEN_type(pg);
        castbuf.nproc = myrank;
        castbuf.psend = pv;
        memcpy(&cast[0],&castbuf,sizeof(coor_comm));           
        memcpy(&cast[sizeof(coor_comm)],msg,sizebuf);          
        free(msg);
        AP_send(buf);
        sendcounts[remoteID]++;    
        delete []tabproc;
      } else {
        int remoteID = (*recup2)[0].first;
        assert(remoteID!=myrank);
        int sizebuf;
        void *msg = de.sendData ((pEntity) pv, remoteID, sizebuf );
        void *buf = AP_alloc(remoteID,de.tag(),sizeof(coor_comm)+sizebuf);
        char *cast = (char *) buf;  
        coor_comm castbuf ; 
        castbuf.X     = P_x(pv);
        castbuf.Y     = P_y(pv);
        castbuf.Z     = P_z(pv);
        pGEntity  pg   = EN_whatIn(pv);
        castbuf.tag   = GEN_tag(pg);
        castbuf.dim   = GEN_type(pg);
        castbuf.nproc = myrank;
        castbuf.psend = pv;
        memcpy(&cast[0],&castbuf,sizeof(coor_comm));           
        memcpy(&cast[sizeof(coor_comm)],msg,sizebuf);          
        free(msg);
        AP_send(buf);
        sendcounts[remoteID]++;    
      }
    }
    VIter_delete(vit);   

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        char * castbuf = (char *) msg;
        coor_comm * coorcom = (coor_comm * )  &castbuf[0];
        pGEntity pent;
        int dim  = coorcom->dim;
        int tag  = coorcom->tag;
        if(dim == 0) {
          pent = (pGEntity) GM_vertexByTag(mesh->model,tag);
        } else if(dim==1) {
          pent = (pGEntity) GM_edgeByTag(mesh->model,tag);     
        } else if(dim==2) {
          pent = (pGEntity) GM_faceByTag(mesh->model,tag);
        } else if(dim==3) {
          pent = (pGEntity) GM_regionByTag(mesh->model,tag);
        } else {
          printf("pbs**** %d\n",dim);
        }
        //         double  X = coorcom->X;
        //         double  Y = coorcom->Y;
        //         double  Z = coorcom->Z;
        //         pVertex pnew = M_createVP(mesh,X,Y,Z,-1,pent);
        double  XYZ[3] = {coorcom->X,coorcom->Y,coorcom->Z};
        pVertex pnew = M_createV2(mesh,XYZ,-1,pent);
        assert(pnew);
        int procdep = coorcom->nproc;
        pVertex pv  = coorcom->psend;
        EN_attachDataInt((pEntity) pnew , tagChange,10);
        VectorOfCommonpVertex_type oldRemoteData;
        oldRemoteData.push_back(CommonpVertex_type(procdep,pv));
        EN_attachDataPtr((pEntity) pnew , tagRemote, 
                         new VectorOfCommonpVertex_type(oldRemoteData));     

        de.receiveData (pnew,from, &castbuf[sizeof(coor_comm)]);
            
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;

    UpdatePointerInterface(mesh);
    UpdatePointerInternal(mesh);
  }

  // -------------------------------------------------------------------
  struct tetra_comm {
    pVertex pdest[4];
    int tag,dim;
  };

  // -------------------------------------------------------------------
  void SendTetra(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de) {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");

    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      int dest;
      int migre = EN_getDataInt((pEntity) pr ,tagElt, &dest);
      if(!migre) continue;
      dest--;
      assert(dest!=myrank);
      pVertex nod[4];
      nod[0] = R_vertex(pr,0);
      nod[1] = R_vertex(pr,1);
      nod[2] = R_vertex(pr,2);
      nod[3] = R_vertex(pr,3);
      int sizebuf;
      void *msg = de.sendData ((pEntity) pr, dest, sizebuf );
      void *buf = AP_alloc(dest,de.tag(),sizeof(tetra_comm)+sizebuf);
      char *cast = (char *) buf;  
      tetra_comm castbuf ; 
      pGEntity pg = EN_whatIn(pr);
      castbuf.tag = GEN_tag(pg);
      castbuf.dim = GEN_type(pg); 
      for(int i = 0 ; i<4 ; i++) {
        int tmp;
        int isChanged = EN_getDataInt((pEntity) nod[i] ,tagChange, &tmp);
        assert(isChanged);
        void *temp_ptr2; 
        int is = EN_getDataPtr((pEntity) nod[i] , tagRemote, &temp_ptr2);
        assert(is);
        VectorOfCommonpVertex_type *recup2 = (VectorOfCommonpVertex_type *) temp_ptr2;
        if((*recup2).size()==1){
          if((tmp==2) || (tmp==10)){
            assert((*recup2)[0].first==dest);
            castbuf.pdest[i] = (*recup2)[0].second;
          } else {
            void *temp_ptr; 
            int is = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
            assert(is);
            const VectorOfCommonpVertex_type *recup = (const VectorOfCommonpVertex_type *) temp_ptr;
            unsigned int j;
            for(j=0 ; j<(*recup).size() ; j++) {
              if((*recup)[j].first == dest) {
                castbuf.pdest[i] = (*recup)[j].second;
                break;
              }
            }
            assert(j<(*recup).size());  
          }
        } else {
          assert((tmp==1) || (tmp==3));
          if(tmp==3) {
            unsigned int j;
            for(j=0 ; j<(*recup2).size() ; j++) {
              if((*recup2)[j].first == dest) {
                castbuf.pdest[i] = (*recup2)[j].second;
                break;
              }
            }
            assert(j<(*recup2).size());
          } else {
            void *temp_ptr; 
            int is = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
            assert(is);
            const VectorOfCommonpVertex_type *recup = (const VectorOfCommonpVertex_type *) temp_ptr;
            unsigned int j;
            for(j=0 ; j<(*recup).size() ; j++) {
              if((*recup)[j].first == dest) {
                castbuf.pdest[i] = (*recup)[j].second;
                break;
              }
            }
            assert(j<(*recup).size());        
          }
        }       
      }  
      memcpy(&cast[0],&castbuf,sizeof(tetra_comm));             
      memcpy(&cast[sizeof(tetra_comm)],msg,sizebuf);
      free(msg);     
      AP_send(buf);
      sendcounts[dest]++;    
    }
    RIter_delete(rit);
  
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        char * castbuf = (char *) msg;
        tetra_comm * regioncom = (tetra_comm * ) &castbuf[0];
        pVertex  p1 = regioncom->pdest[0];
        int tmp;
        int isChanged = EN_getDataInt((pEntity) p1 ,tagChange, &tmp);
        assert(isChanged);      
        pVertex  p2 = regioncom->pdest[1];
        isChanged = EN_getDataInt((pEntity) p2 ,tagChange, &tmp);
        assert(isChanged);
        pVertex  p3 = regioncom->pdest[2];
        isChanged = EN_getDataInt((pEntity) p3 ,tagChange, &tmp);
        assert(isChanged);
        pVertex  p4 = regioncom->pdest[3];
        isChanged = EN_getDataInt((pEntity) p4 ,tagChange, &tmp);
        assert(isChanged);
        pFace    pface[4];
        pGEntity pg;
        int dim  = regioncom->dim;
        int tag  = regioncom->tag;
        if(dim==2) {
          pg = (pGEntity) GM_faceByTag(mesh->model,tag);
        } else if(dim==3) {
          pg = (pGEntity) GM_regionByTag(mesh->model,tag);
        } else {
          printf("----pbs faces**** %d\n",dim);
        }
        //         pface[0] =  F_exist(2,p1,p2,p3,0);
        //         assert(pface[0]);
        //         pface[1] =  F_exist(2,p1,p2,p4,0);
        //         assert(pface[1]);
        //         pface[2] =  F_exist(2,p2,p3,p4,0);
        //         assert(pface[2]);
        //         pface[3] =  F_exist(2,p1,p3,p4,0);
        //         assert(pface[3]);
        
        pface[0] =  F_exist(p1,p2,p3,0);
        assert(pface[0]);
        pface[1] =  F_exist(p1,p2,p4,0);
        assert(pface[1]);
        pface[2] =  F_exist(p2,p3,p4,0);
        assert(pface[2]);
        pface[3] =  F_exist(p1,p3,p4,0);
        assert(pface[3]);

        pRegion pr = M_createR(mesh,4,pface,pg);
        assert(pr);
        de.receiveData (pr,from, &castbuf[sizeof(tetra_comm)]);
             
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;

  }
  // -------------------------------------------------------------------
  struct face_comm2 {
    pVertex pdest[3];
    int tag,dim;
  };

  // -------------------------------------------------------------------
  void SendFaces(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de) {
    int Dim = (mesh->tets.empty()) ? 2 : 3;
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");

    if(Dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        int dest;
        int migre = EN_getDataInt((pEntity) pface ,tagElt, &dest);
        if(!migre) continue;
        dest--;
        assert(dest!=myrank);
        pVertex nod[3];
        pface->getNodes(nod);
        int sizebuf;
        void *msg = de.sendData ((pEntity) pface, dest, sizebuf );
        void *buf = AP_alloc(dest,de.tag(),sizeof(face_comm2) +sizebuf);
        char *cast = (char *) buf;  
       
        face_comm2 castbuf; 
        pGEntity pg = EN_whatIn(pface);
        castbuf.tag = GEN_tag(pg);
        castbuf.dim = GEN_type(pg); 
        for(int i = 0 ; i<3 ; i++) {
          int tmp;
          int isChanged = EN_getDataInt((pEntity) nod[i] ,tagChange, &tmp);
          assert(isChanged);
          void *temp_ptr2; 
          int is = EN_getDataPtr((pEntity) nod[i] , tagRemote, &temp_ptr2);
          assert(is);
          VectorOfCommonpVertex_type *recup2 = 
            (VectorOfCommonpVertex_type *) temp_ptr2;
          if((*recup2).size()==1){
            if((tmp==2) || (tmp==10)){
              assert((*recup2)[0].first==dest);
              castbuf.pdest[i] = (*recup2)[0].second;
            } else {
              void *temp_ptr; 
              int is = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
              assert(is);
              const VectorOfCommonpVertex_type *recup = 
                (const VectorOfCommonpVertex_type *) temp_ptr;
              unsigned int j;
              for(j=0 ; j<(*recup).size() ; j++) {
                if((*recup)[j].first == dest) {
                  castbuf.pdest[i] = (*recup)[j].second;
                  break;
                }
              }
              assert(j<(*recup).size());    
            }
          } else {
            assert((tmp==1) || (tmp==3));
            if(tmp==3) {
              unsigned int j;
              for(j=0 ; j<(*recup2).size() ; j++) {
                if((*recup2)[j].first == dest) {
                  castbuf.pdest[i] = (*recup2)[j].second;
                  break;
                }
              }
              assert(j<(*recup2).size());
            } else {
              void *temp_ptr; 
              int is = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
              assert(is);
              const VectorOfCommonpVertex_type *recup = 
                (const VectorOfCommonpVertex_type *) temp_ptr;
              unsigned int j;
              for(j=0 ; j<(*recup).size() ; j++) {
                if((*recup)[j].first == dest) {
                  castbuf.pdest[i] = (*recup)[j].second;
                  break;
                }
              }
              assert(j<(*recup).size());          
            }
          } 
        }    
        memcpy(&cast[0],&castbuf,sizeof(face_comm2));         
        memcpy(&cast[sizeof(face_comm2)],msg,sizebuf);
        free(msg);      
        AP_send(buf);
        sendcounts[dest]++;    
      }
      FIter_delete(fit);
    } else {//Dim==3
      if(Dim!=3) throw;
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        pPList prlist = F_regions(pface);
        pRegion pr;
        void* iter=0;
        while( (pr =(pRegion) PList_next(prlist,&iter))){
          int dest;
          int migre = EN_getDataInt((pEntity) pr ,tagElt, &dest);
          if(!migre) continue;
          dest--;
          assert(dest!=myrank);
          ((MDB_Triangle*)pface)->del((MDB_Tet*)pr);   

          pVertex nod[3];
          pface->getNodes(nod);
          int sizebuf;
          void *msg = de.sendData ((pEntity) pface, dest, sizebuf );
          void *buf = AP_alloc(dest,de.tag(),sizeof(face_comm2)+sizebuf);
          char *cast = (char*) buf;
          face_comm2 castbuf; 
          pGEntity pg = EN_whatIn(pface);
          castbuf.tag = GEN_tag(pg);
          castbuf.dim = GEN_type(pg); 
          for(int i = 0 ; i<3 ; i++) {
            int tmp;
            int isChanged = EN_getDataInt((pEntity) nod[i] ,tagChange, &tmp);
            assert(isChanged);
            void *temp_ptr2; 
            int is = EN_getDataPtr((pEntity) nod[i] , tagRemote, &temp_ptr2);
            assert(is);
            VectorOfCommonpVertex_type *recup2 = 
              (VectorOfCommonpVertex_type *) temp_ptr2;
            if((*recup2).size()==1){
              if((tmp==2) || (tmp==10)){
                assert((*recup2)[0].first==dest);
                castbuf.pdest[i] = (*recup2)[0].second;
              } else {
                void *temp_ptr; 
                int is = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
                assert(is);
                const VectorOfCommonpVertex_type *recup = 
                  (const VectorOfCommonpVertex_type *) temp_ptr;
                unsigned int j;
                for(j=0 ; j<(*recup).size() ; j++) {
                  if((*recup)[j].first == dest) {
                    castbuf.pdest[i] = (*recup)[j].second;
                    break;
                  }
                }
                assert(j<(*recup).size());    
              }
            } else {
              assert((tmp==1) || (tmp==3));
              if(tmp==3) {
                unsigned int j;
                for(j=0 ; j<(*recup2).size() ; j++) {
                  if((*recup2)[j].first == dest) {
                    castbuf.pdest[i] = (*recup2)[j].second;
                    break;
                  }
                }
                assert(j<(*recup2).size());
              } else {
                void *temp_ptr; 
                int is = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
                assert(is);
                const VectorOfCommonpVertex_type *recup = 
                  (const VectorOfCommonpVertex_type *) temp_ptr;
                unsigned int j;
                for(j=0 ; j<(*recup).size() ; j++) {
                  if((*recup)[j].first == dest) {
                    castbuf.pdest[i] = (*recup)[j].second;
                    break;
                  }
                }
                assert(j<(*recup).size());          
              }
            } 
          }    
          memcpy(&cast[0],&castbuf,sizeof(face_comm2));           
          memcpy(&cast[sizeof(face_comm2)],msg,sizebuf);
          free(msg);  
          AP_send(buf);
          sendcounts[dest]++;    
        }
      }
      FIter_delete(fit);  
    }
  
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        char * castbuf = (char *) msg;

        face_comm2 * facecom = (face_comm2 * ) &castbuf[0];
        pVertex  p1 = facecom->pdest[0];
        int tmp;
        int isChanged = EN_getDataInt((pEntity) p1 ,tagChange, &tmp);
        assert(isChanged);      
        pVertex  p2 = facecom->pdest[1];
        isChanged = EN_getDataInt((pEntity) p2 ,tagChange, &tmp);
        assert(isChanged);
        pVertex  p3 = facecom->pdest[2];
        isChanged = EN_getDataInt((pEntity) p3 ,tagChange, &tmp);
        assert(isChanged);
        pEdge    pe[3];
        pGEntity pg;
        int dim  = facecom->dim;
        int tag  = facecom->tag;
        if(dim==2) {
          pg = (pGEntity) GM_faceByTag(mesh->model,tag);
        } else if(dim==3) {
          pg = (pGEntity) GM_regionByTag(mesh->model,tag);
        } else {
          printf("----pbs faces**** %d\n",dim);
        }
        pe[0] =  E_exist(p1,p2);
        assert(pe[0]);
        pe[1] =  E_exist(p2,p3);
        assert(pe[1]);
        pe[2] =  E_exist(p1,p3);
        assert(pe[2]);
       
        //      pFace pface =F_exist(2,pe[0],pe[1],pe[2],0);
        //      if(!pface) {
        pFace  pface = M_createF(mesh,3,pe,pg);
        assert(pface);
        de.receiveData (pface,from, &castbuf[sizeof(face_comm2)]);
       
        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;

  }

  // -------------------------------------------------------------------
  struct edge_comm {
    pVertex pdest[2];
    int tag,dim;
  };

  // -------------------------------------------------------------------
  void SendEdges(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de) {
    int Dim = (mesh->tets.empty()) ? 2 : 3;
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
  
    if(Dim==2) {
      EIter eit = M_edgeIter(mesh);
      pEdge pedge;  
      while ((pedge = EIter_next(eit))) {
        pVertex p[2];
        p[0] = pedge->p1;
        p[1] = pedge->p2;
        int tmp1;
        int isChanged1 = EN_getDataInt((pEntity) p[0] ,tagChange, &tmp1);
        int tmp2;
        int isChanged2 = EN_getDataInt((pEntity) p[1] ,tagChange, &tmp2);
        if(!isChanged1 || !isChanged2) continue;
        int nbfaces = pedge->numfaces();
        for(int k = 0 ; k < nbfaces ; k++) {
          pFace pface = pedge->faces(k);
          int dest;
          int migre = EN_getDataInt((pEntity) pface ,tagElt, &dest);
          if(!migre) continue;
          dest--;
          assert(dest!=myrank);
          pedge->del((MDB_Triangle*)pface);  
          int sizebuf;
          void *msg = de.sendData ((pEntity) pedge, dest, sizebuf );
    
          void *buf = AP_alloc(dest,de.tag(),sizeof(edge_comm) + sizebuf);
          char *cast = (char *) buf;  
          edge_comm castbuf; 
          pGEntity  pg = EN_whatIn(pedge);
          castbuf.tag = GEN_tag(pg);
          castbuf.dim = GEN_type(pg);
          for(int i = 0 ; i< 2 ; i++) {
            int tmp;
            int isChanged = EN_getDataInt((pEntity) p[i] ,tagChange, &tmp);
            assert(isChanged);
            void *temp_ptr2; 
            int is = EN_getDataPtr((pEntity) p[i] , tagRemote, &temp_ptr2);
            assert(is);
            VectorOfCommonpVertex_type *recup2 = 
              (VectorOfCommonpVertex_type *) temp_ptr2;
            if((*recup2).size()==1){
              if((tmp==2) || (tmp==10)){
                assert((*recup2)[0].first==dest);
                castbuf.pdest[i] = (*recup2)[0].second;
              } else {
                void *temp_ptr; 
                int is = EN_getDataPtr((pEntity) p[i] , tagData, &temp_ptr);
                assert(is);
                const VectorOfCommonpVertex_type *recup = 
                  (const VectorOfCommonpVertex_type *) temp_ptr;
                unsigned int j;
                for(j=0 ; j<(*recup).size() ; j++) {
                  if((*recup)[j].first == dest) {
                    castbuf.pdest[i] = (*recup)[j].second;
                    break;
                  }
                }
                assert(j<(*recup).size());  
              }
            } else {
              assert((tmp==1) || (tmp==3));
              if(tmp==3) {
                unsigned int j;
                for(j=0 ; j<(*recup2).size() ; j++) {
                  if((*recup2)[j].first == dest) {
                    castbuf.pdest[i] = (*recup2)[j].second;
                    break;
                  }
                }
                assert(j<(*recup2).size());
              } else {
                void *temp_ptr; 
                int is = EN_getDataPtr((pEntity) p[i] , tagData, &temp_ptr);
                assert(is);
                const VectorOfCommonpVertex_type *recup = 
                  (const VectorOfCommonpVertex_type *) temp_ptr;
                unsigned int j;
                for(j=0 ; j<(*recup).size() ; j++) {
                  if((*recup)[j].first == dest) {
                    castbuf.pdest[i] = (*recup)[j].second;
                    break;
                  }
                }
                assert(j<(*recup).size());        
              }
            }   
          }
          memcpy(&cast[0],&castbuf,sizeof(edge_comm));           
          memcpy(&cast[sizeof(edge_comm)],msg,sizebuf);
          free(msg);  
          AP_send(buf);
          sendcounts[dest]++;      
        }
      } 
      EIter_delete(eit); 
    } else {//Dim==3
      if(Dim!=3) throw;
      EIter eit = M_edgeIter(mesh);
      pEdge pedge;  
      while ((pedge = EIter_next(eit))) {
        pVertex p[2];
        p[0] = pedge->p1;
        p[1] = pedge->p2;
        int tmp1;
        int isChanged1 = EN_getDataInt((pEntity) p[0] ,tagChange, &tmp1);
        int tmp2;
        int isChanged2 = EN_getDataInt((pEntity) p[1] ,tagChange, &tmp2);
        if(!isChanged1 || !isChanged2) continue;
        pPList prlist = E_regions(pedge);
        pRegion pr;
        void* iter=0;
        while( (pr =(pRegion) PList_next(prlist,&iter))){
          int dest;
          int migre = EN_getDataInt((pEntity) pr ,tagElt, &dest);
          if(!migre) continue;
          dest--;
          assert(dest!=myrank);  
          int sizebuf;
          void *msg = de.sendData ((pEntity) pedge, dest, sizebuf );
      
          void *buf = AP_alloc(dest,de.tag(),sizeof(edge_comm) +sizebuf);
          char *cast = (char *) buf;  
 
          edge_comm castbuf; 
          pGEntity  pg = EN_whatIn(pedge);
          castbuf.tag = GEN_tag(pg);
          castbuf.dim = GEN_type(pg);       
        
          for(int i = 0 ; i< 2 ; i++) {
            int tmp;
            int isChanged = EN_getDataInt((pEntity) p[i] ,tagChange, &tmp);
            assert(isChanged);
            void *temp_ptr2; 
            int is = EN_getDataPtr((pEntity) p[i] , tagRemote, &temp_ptr2);
            assert(is);
            VectorOfCommonpVertex_type *recup2 = 
              (VectorOfCommonpVertex_type *) temp_ptr2;
            if((*recup2).size()==1){
              if((tmp==2) || (tmp==10)){
                assert((*recup2)[0].first==dest);
                castbuf.pdest[i] = (*recup2)[0].second;
              } else {
                void *temp_ptr; 
                int is = EN_getDataPtr((pEntity) p[i] , tagData, &temp_ptr);
                assert(is);
                const VectorOfCommonpVertex_type *recup = 
                  (const VectorOfCommonpVertex_type *) temp_ptr;
                unsigned int j;
                for(j=0 ; j<(*recup).size() ; j++) {
                  if((*recup)[j].first == dest) {
                    castbuf.pdest[i] = (*recup)[j].second;
                    break;
                  }
                }
                assert(j<(*recup).size());  
              }
            } else {
              assert((tmp==1) || (tmp==3));
              if(tmp==3) {
                unsigned int j;
                for(j=0 ; j<(*recup2).size() ; j++) {
                  if((*recup2)[j].first == dest) {
                    castbuf.pdest[i] = (*recup2)[j].second;
                    break;
                  }
                }
                assert(j<(*recup2).size());
              } else {
                void *temp_ptr; 
                int is = EN_getDataPtr((pEntity) p[i] , tagData, &temp_ptr);
                assert(is);
                const VectorOfCommonpVertex_type *recup = 
                  (const VectorOfCommonpVertex_type *) temp_ptr;
                unsigned int j;
                for(j=0 ; j<(*recup).size() ; j++) {
                  if((*recup)[j].first == dest) {
                    castbuf.pdest[i] = (*recup)[j].second;
                    break;
                  }
                }
                assert(j<(*recup).size());        
              }
            }   
          } 
          memcpy(&cast[0],&castbuf,sizeof(edge_comm));           
          memcpy(&cast[sizeof(edge_comm)],msg,sizebuf);
          free(msg);  
          AP_send(buf);
          sendcounts[dest]++;  
        } 
      }   
      EIter_delete(eit);
    }
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);

    /*receive pointers*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        char * castbuf = (char *) msg;
        edge_comm * edgecom = (edge_comm * ) &castbuf[0];
        pVertex  p1 = edgecom->pdest[0];
        int tmp;
        int isChanged = EN_getDataInt((pEntity) p1 ,tagChange, &tmp);
        assert(isChanged);      
        pVertex  p2 = edgecom->pdest[1];
        isChanged = EN_getDataInt((pEntity) p2 ,tagChange, &tmp);
        assert(isChanged);
        pGEntity pg;
        int dim  = edgecom->dim;
        int tag  = edgecom->tag;
        if(dim==1) {
          pg = (pGEntity) GM_edgeByTag(mesh->model,tag);
        } else if(dim==2) {
          pg = (pGEntity) GM_faceByTag(mesh->model,tag);
        } else if(dim==3) {
          pg = (pGEntity) GM_regionByTag(mesh->model,tag);
        } else {
          printf("----pbs**** %d\n",dim);
        }
        pEdge     e = E_exist(p1,p2);
        if (!e)   e = M_createE(mesh,p1,p2,pg);
        de.receiveData (e,from, &castbuf[sizeof(edge_comm)]);

        AP_free(msg);      
      }
    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;


  }

  // -------------------------------------------------------------------
  void SendElt(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de) {

    SendEdges(mesh,tagElt,de);
  
    SendFaces(mesh,tagElt,de);
  
    int dim = (mesh->tets.empty()) ? 2 : 3;
    if(dim==3) SendTetra(mesh,tagElt,de);
  }

  // -------------------------------------------------------------------
  void DeleteEntities(pMesh mesh,pMeshDataId tagElt) {
    int nproc,myrank;
    int Dim = (mesh->tets.empty()) ? 2 : 3;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);   

    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
  
    if(Dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        int dest;
        int migre = EN_getDataInt((pEntity) pface ,tagElt, &dest);
        if(!migre) continue;
        dest--;
        assert(dest!=myrank);
        EN_deleteData((pEntity)pface, tagElt);
        M_removeFace(mesh,pface);
      }
      FIter_delete(fit);
      MD_deleteMeshDataId(tagElt);
    } else {
      RIter rit = M_regionIter(mesh);
      pRegion pr;  
      while ((pr = RIter_next(rit))) {
        int dest;
        int migre = EN_getDataInt((pEntity) pr ,tagElt, &dest);
        if(!migre) continue;
        dest--;
        assert(dest!=myrank);
        EN_deleteData((pEntity)pr, tagElt);
        M_removeRegion(mesh,pr);
      }
      RIter_delete(rit);    
      MD_deleteMeshDataId(tagElt);  
    
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        int num = F_numRegions(pface);
        if(!num) M_removeFace(mesh,pface);
      }
      FIter_delete(fit);
    }
  
    EIter eit = M_edgeIter(mesh);
    pEdge ped;  
    while ((ped = EIter_next(eit))) {
      int num = E_numFaces(ped);
      if(!num)  M_removeEdge(mesh,ped);
    }
    EIter_delete(eit);

    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      int tmp;
      int isChanged = EN_getDataInt((pEntity) pv ,tagChange, &tmp);
      if(!isChanged) continue;
      void *temp_ptr2; 
      int is2 = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr2);
#ifdef DEBUG
      if(tmp==10)assert(!is2);
      if(tmp==1) assert(is2);
#endif
      void *temp_ptr; 
      int is = EN_getDataPtr((pEntity) pv , tagRemote, &temp_ptr);
    
      
      int num = V_numEdges(pv);
      if(tmp==3 ) assert(!num); 
      if(tmp==11) {
        assert(num);
        assert(is2);
        EN_deleteData((pEntity)pv, tagData);
      }  
      if(tmp==2)  assert(!num); 
      if(tmp==1) assert(num);
      if(!num) assert(!V_numFaces(pv));
    
      if(is) EN_deleteData((pEntity)pv,tagRemote);
      int isMaster = EN_getDataInt((pEntity) pv ,tagMaster, &tmp);

      if(isMaster)  EN_deleteData((pEntity)pv, tagMaster);
      EN_deleteData((pEntity)pv, tagChange);

      if(!num) M_removeVertex(mesh,pv);
    }
    VIter_delete(vit);
    MD_deleteMeshDataId(tagChange);
    MD_deleteMeshDataId(tagRemote);  
    MD_deleteMeshDataId(tagMaster);  
  }
  // -------------------------------------------------------------------
  void DeleteEntitiesAndData(pMesh mesh, pMeshDataId tagElt, 
                             MDB_DataExchanger &de )
  {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);   
    pMeshDataId tagRemote = MD_lookupMeshDataId("RemoteStructure");
    pMeshDataId tagMaster = MD_lookupMeshDataId("isMaster");
    pMeshDataId tagChange = MD_lookupMeshDataId("IsChange");
    int Dim = (mesh->tets.empty()) ? 2 : 3; 
    if(Dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        EN_deleteData((pEntity)pface, tagRemote);
        EN_deleteData((pEntity)pface, tagMaster);
        EN_deleteData((pEntity)pface, tagChange);
        int dest;
        int migre = EN_getDataInt((pEntity) pface ,tagElt, &dest);
        if(!migre) continue;
        dest--;
        assert(dest!=myrank);
        EN_deleteData((pEntity)pface, tagElt);
        de.deleteExternalData( (pEntity) pface );
        M_removeFace(mesh,pface);
      }
      FIter_delete(fit);
      MD_deleteMeshDataId(tagElt);
    } 
    else {
      RIter rit = M_regionIter(mesh);
      pRegion pr;  
      while ((pr = RIter_next(rit))) {
        EN_deleteData((pEntity)pr, tagRemote);
        EN_deleteData((pEntity)pr, tagMaster);
        EN_deleteData((pEntity)pr, tagChange);
        int dest;
        int migre = EN_getDataInt((pEntity) pr ,tagElt, &dest);
        if(!migre) continue;
        dest--;
        assert(dest!=myrank);
        EN_deleteData((pEntity)pr, tagElt);
        de.deleteExternalData( (pEntity) pr );
        M_removeRegion(mesh,pr);
      }
      RIter_delete(rit);    
      MD_deleteMeshDataId(tagElt);  
    
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
      
        EN_deleteData((pEntity)pface, tagRemote);
        EN_deleteData((pEntity)pface, tagMaster);
        EN_deleteData((pEntity)pface, tagChange);
        int num = F_numRegions(pface);
        if(!num){
          de.deleteExternalData( (pEntity) pface );
          M_removeFace(mesh,pface);
        }
      }
      FIter_delete(fit);
    }
  
    EIter eit = M_edgeIter(mesh);
    pEdge ped;  
    while ((ped = EIter_next(eit))) {
      EN_deleteData((pEntity)ped, tagRemote);
      EN_deleteData((pEntity)ped, tagMaster);
      EN_deleteData((pEntity)ped, tagChange);
      int num = E_numFaces(ped);
      if(!num){
        de.deleteExternalData( (pEntity) ped );      
        M_removeEdge(mesh,ped);
      }
    }
    EIter_delete(eit);
  
    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      EN_deleteData((pEntity)pv, tagRemote);
      EN_deleteData((pEntity)pv, tagMaster);
      EN_deleteData((pEntity)pv, tagChange);
      int num = V_numEdges(pv);
      if(!num){
        de.deleteExternalData( (pEntity) pv );   
        M_removeVertex(mesh,pv);   
      }
    }
    VIter_delete(vit);
    MD_deleteMeshDataId(tagChange);
    MD_deleteMeshDataId(tagRemote);  
    MD_deleteMeshDataId(tagMaster);  
  
  }
  // -------------------------------------------------------------------
  void loadBalancing2(pMesh mesh, pMeshDataId tagElt, MDB_DataExchanger &de)
  {
    // 0 the destination of local subentitiees 
    MarkEltSubEntities(mesh, tagElt);
    // 1 entities Migration
    pMeshDataId tagDest = MD_lookupMeshDataId("tagDestinations");
    MigrateEntitiesAndData(mesh, tagDest, de);
    // 2 delete entities
    DeleteEntitiesAndData( mesh,tagElt, de );
 
    pMeshDataId tagVertex = MD_lookupMeshDataId("RemotePoint");
    MD_deleteMeshDataId(tagVertex);
    MD_deleteMeshDataId(tagDest);
    mesh->initializeIdData();
  }
  // -------------------------------------------------------------------
  void loadBalancing(pMesh mesh, pMeshDataId tagElt, MDB_DataExchanger &de)
  {
#ifdef DEBUG
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    checkRemotePointer(mesh,tagData);
    puts("check first ok");
#endif
    // 1) Mark local and distant vertices with their new proc(s)
    MarkEltVertex(mesh,tagElt);

#ifdef DEBUG  
    checkRemotePointer(mesh,tagData);
    puts("check second ok");
#endif
  
    // 2) update and send new interfaces points
    UpdateInterfaces(mesh,de);

#ifdef DEBUG
    checkRemotePointer2(mesh,tagData);
    puts("check fifth ok");
#endif
  
    // 3) send other points and update RemotePointer
    SendVertex(mesh,de);
  
    // 4) send elt
    SendElt(mesh,tagElt,de);
  
    //  // 5) delete entities
    //  DeleteEntities(mesh,tagElt);
    DeleteEntitiesAndData( mesh,tagElt, de );

    //  // 6) classification 
    mesh->classify_unclassified_entities();

 
    pMeshDataId tagVertex = MD_lookupMeshDataId("RemotePoint");
    MD_deleteMeshDataId(tagVertex);
#ifdef PARALLEL

    // ----------------------------------------------
    // ------ Tagging inter-partition nodes
    // ----------------------------------------------
  
    pMeshDataId tag = MD_lookupMeshDataId("RemotePoint");
  
    V_createInfoInterface(mesh,tag);
    E_createInfoInterface(mesh,tag);
    F_createInfoInterface(mesh,tag);
    mesh->initializeIdData();
  
#endif

#ifdef DEBUG  
    checkRemotePointer(mesh,tagData);
    puts("last check ok");
#endif  

  }
  // -------------------------------------------------------------------

} // End of namespace MAd

#endif
