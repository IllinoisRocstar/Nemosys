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
#include <cassert>

#include "mpi.h"
#include "autopack.h"

#include "MeshDataBase.h"
#include "MeshDataBaseParallelInterface.h"
#include "CheckOrientation.h"

namespace MAd {
  
  struct edge_comm
  {
    pVertex p1,p2;                //pointer in dest
    int     sendID;       
  };
  
  struct face_comm {
    pVertex pdest[3];
  };  

  int CheckEdgesOrientation(pMesh mesh){
#ifdef DEBUG
    int norient = 0;
#endif
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;

    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    EIter eit = M_edgeIter(mesh);
    pEdge pedge;  
    while ((pedge = EIter_next(eit))) {
      pVertex p[2];
      p[0] = E_vertex(pedge,0);
      p[1] = E_vertex(pedge,1);
      int tmp1;
      int isInterface1 = EN_getDataInt((pEntity) p[0] ,tagData, &tmp1);
      int tmp2;
      int isInterface2 = EN_getDataInt((pEntity) p[1] ,tagData, &tmp2);
      if(!isInterface1 || !isInterface2) continue;
      const std::vector<std::pair<int , pVertex> > *recup1 = (std::vector<std::pair<int , pVertex> > *) tmp1;
      const std::vector<std::pair<int , pVertex> > *recup2 = (std::vector<std::pair<int , pVertex> > *) tmp2;
      int size1 = (*recup1).size();
      int size2 = (*recup2).size();
      assert(size1);assert(size2);
      int *tab=new int[size1];
      for(int i=0 ; i< size1 ; i++) tab[i] = (*recup1)[i].first;
      //check if myrank must send
      int nSender = myrank;
      for(int j=0 ; j< size2 ; j++) {    
        int iProc = (*recup2)[j].first;
        int i;
        for(i=0 ; i<size1 ; i++) {
          if(iProc == tab[i]) break;
        }
        if(i<size1) {
          if(iProc < nSender) nSender = iProc;
        }
      }
      if(nSender != myrank) continue;
      for(int j=0 ; j< size2 ; j++) {    
        int iProc = (*recup2)[j].first;
        int i;
        for(i=0 ; i<size1 ; i++) {
          if(iProc == tab[i]) break;
        }
        if(i < size1) {  
          pVertex remote1 = (*recup1)[i].second;
          pVertex remote2 = (*recup2)[j].second;
          assert(iProc != myrank);
          assert(iProc > myrank);
          void *buf = AP_alloc(iProc,444,sizeof(edge_comm));
          edge_comm *castbuf = (edge_comm *) buf;
          castbuf->p1	   = remote1;
          castbuf->p2	   = remote2,
            castbuf->sendID    = myrank;
          AP_send(buf);
          sendcounts[iProc]++;         
        }
      }  
      delete []tab;
    }
    EIter_delete(eit); 
    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    int message=0;
    int count;
  
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
        edge_comm * castbuf = (edge_comm*) msg;
        pVertex p1recv,p2recv;
        p1recv = castbuf -> p1;
        p2recv = castbuf -> p2;
        assert(p1recv);assert(p2recv);
        int nprocrecv = castbuf->sendID;
        assert(nprocrecv==from);
        pEdge pe = E_exist(p1recv,p2recv);
        if(pe) {
          //check orientation
          if(E_vertex(pe,0) != p1recv) {
#ifdef DEBUG
            norient++;
#endif
            pe->p1 = p1recv;
            pe->p2 = p2recv;
            assert(E_vertex(pe,0) == p1recv);
          }
        } 
        AP_free(msg);
      }
    }    
    AP_check_sends(AP_WAITALL); 
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
#ifdef DEBUG
    return norient;
#else    
    return 0;
#endif  
  }


  int CheckFacesOrientation(pMesh mesh){
#ifdef DEBUG
    int norient = 0;
#endif
    int mysize,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
 
    int *sendcounts = new int[mysize];
    for(int i=0;i<mysize;i++) sendcounts[i]=0;
  
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
  
    FIter fit = M_faceIter(mesh);
    pFace pface;  
    while ((pface = FIter_next(fit))) {
      pVertex p1 = F_vertex(pface,0);
      pVertex p2 = F_vertex(pface,1);
      pVertex p3 = F_vertex(pface,2);
      void *temp_ptr1,*temp_ptr2,*temp_ptr3; 
      int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
      int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
      int isInterface3 = EN_getDataPtr((pEntity) p3 , tagData, &temp_ptr3);
      if(!(isInterface1 && isInterface2 && isInterface3)) continue;
   
      const std::vector<std::pair<int , pVertex> > *recup1 = (std::vector<std::pair<int , pVertex> > *) temp_ptr1;
      const std::vector<std::pair<int , pVertex> > *recup2 = (std::vector<std::pair<int , pVertex> > *) temp_ptr2;
      const std::vector<std::pair<int , pVertex> > *recup3 = (std::vector<std::pair<int , pVertex> > *) temp_ptr3;
      int size = (*recup1).size();
      int *tab=new int[size];
      for(int i=0 ; i< size ; i++) tab[i] = (*recup1)[i].first;
    
      int size2 = (*recup2).size();
      int *tab2=new int[size2];
      for(int i=0 ; i< size2 ; i++) tab2[i] = (*recup2)[i].first;

      // check if I must send
      int nSender = myrank;
      for(unsigned int k=0 ; k<(*recup3).size() ; k++) {	 
        int iProc = (*recup3)[k].first;
        int i;
        for(i=0 ; i<size ; i++) {
          if(iProc == tab[i]) break;
        }
        int j;
        for(j=0 ; j<size2 ; j++) {
          if(iProc == tab2[j]) break;
        }       
        if(i < size && j < size2) {
          if(iProc < nSender) nSender = iProc;
        }
      }
    
      if(nSender != myrank) continue;
    
      for(unsigned int k=0 ; k<(*recup3).size() ; k++) {	 
        int iProc = (*recup3)[k].first;
        int i;
        for(i=0 ; i<size ; i++) {
          if(iProc == tab[i]) break;
        }
        int j;
        for(j=0 ; j<size2 ; j++) {
          if(iProc == tab2[j]) break;
        }       
        if(i < size && j < size2) {
          assert(tab[i]==iProc);
          assert(tab2[j]==iProc);
        
          pVertex remote1 = (*recup1)[i].second;
          pVertex remote2 = (*recup2)[j].second;
          pVertex remote3 = (*recup3)[k].second;
        
          void *buf = AP_alloc(iProc,444,sizeof(face_comm));       
          face_comm *castbuf = (face_comm *) buf;
          castbuf->pdest[0]	   = remote1;
          castbuf->pdest[1]	   = remote2,
            castbuf->pdest[2]	   = remote3,
            AP_send(buf);
          sendcounts[iProc]++;    
        }
      }
      delete []tab;
      delete []tab2;
    }      
    FIter_delete(fit);  
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
        face_comm * facecom = (face_comm * ) msg;
        pVertex  p1 = facecom->pdest[0];
        pVertex  p2 = facecom->pdest[1];
        pVertex  p3 = facecom->pdest[2];
        pEdge    pe[3];
        pe[0] =  E_exist(p1,p2);
        pe[1] =  E_exist(p2,p3);
        pe[2] =  E_exist(p1,p3);
      
        if(pe[0] && pe[1] && pe[2]) {
          // pFace pface = F_exist(2,p1,p2,p3,0);
          pFace pface = F_exist(p1,p2,p3,0);
          assert(pface);
      
          //check orientation
          if(F_edge(pface,0) != pe[0] || F_edge(pface,1) != pe[1]) {
            ((MDB_Triangle*)pface)->e1 = pe[0];
            ((MDB_Triangle*)pface)->e2 = pe[1];
            ((MDB_Triangle*)pface)->e3 = pe[2];
#ifdef DEBUG
            norient++;
#endif	  	
          }        
        }       
        AP_free(msg);      
      }    
    }      
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
#ifdef DEBUG
    return norient;
#else  
    return 0;
#endif	  	
  }

/*
// version with sets for the parallel nodes
int CheckEdgesOrientation(pMesh mesh){
#ifdef DEBUG
int norient = 0;
#endif
int nproc,myrank;
MPI_Comm_size(MPI_COMM_WORLD, &nproc);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

int *sendcounts = new int[nproc];
for(int i=0;i<nproc;i++)sendcounts[i]=0;

pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
EIter eit = M_edgeIter(mesh);
pEdge pedge;  
while ((pedge = EIter_next(eit))) {
pVertex p[2];
p[0] = E_vertex(pedge,0);
p[1] = E_vertex(pedge,1);
int tmp1;
int isInterface1 = EN_getDataInt((pEntity) p[0] ,tagData, &tmp1);
int tmp2;
int isInterface2 = EN_getDataInt((pEntity) p[1] ,tagData, &tmp2);
if(!isInterface1 || !isInterface2) continue;
const std::set<std::pair<int, MDB_Point*> > *recup1 = (std::set<std::pair<int, MDB_Point*> > *) temp_ptr1;
const std::set<std::pair<int, MDB_Point*> > *recup2 = (std::set<std::pair<int, MDB_Point*> > *) temp_ptr2;
int size1 = (*recup1).size();
int size2 = (*recup2).size();
assert(size1);assert(size2);
int *tab=new int[size1];
std::set<std::pair<int, MDB_Point*> >::const_iterator rIter1 = recup1->begin();
for (int i=0; rIter1 != recup1->end(); rIter1++, i++) tab[i] = (*rIter1).first;
int nSender = myrank;
std::set<std::pair<int, MDB_Point*> >::const_iterator rIter2 = recup2->begin();
for (int j=0; rIter2 != recup2->end(); rIter2++, j++) {
int iProc = (*rIter2).first;
int i;
for(i=0 ; i<size1 ; i++) {
if(iProc == tab[i]) break;
}
if(i<size1) {
if(iProc < nSender) nSender = iProc;
}
}
if(nSender != myrank) continue;
rIter2 = recup2->begin();
for (int j=0; rIter2 != recup2->end(); rIter2++, j++) {
int iProc = (*rIter2).first;
int i;
rIter1 = recup1->begin();
for (i=0; rIter1 != recup1->end(), i<size1; rIter1++, i++) {
if(iProc == tab[i]) break;
}
if(i < size1) {
pVertex remote1 = (*rIter1).second;
pVertex remote2 = (*rIter2).second;
assert(iProc != myrank);
assert(iProc > myrank);
void *buf = AP_alloc(iProc,444,sizeof(edge_comm));
edge_comm *castbuf = (edge_comm *) buf;
castbuf->p1	   = remote1;
castbuf->p2	   = remote2,
castbuf->sendID    = myrank;
AP_send(buf);
sendcounts[iProc]++;         
}
}  
delete []tab;
}
EIter_delete(eit); 
AP_check_sends(AP_NOFLAGS);
AP_reduce_nsends(sendcounts);
  
int message=0;
int count;
  
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
edge_comm * castbuf = (edge_comm*) msg;
pVertex p1recv,p2recv;
p1recv = castbuf -> p1;
p2recv = castbuf -> p2;
assert(p1recv);assert(p2recv);
int nprocrecv = castbuf->sendID;
assert(nprocrecv==from);
pEdge pe = E_exist(p1recv,p2recv);
if(pe) {
//check orientation
if(E_vertex(pe,0) != p1recv) {
#ifdef DEBUG
norient++;
#endif
pe->p1 = p1recv;
pe->p2 = p2recv;
assert(E_vertex(pe,0) == p1recv);
}
} 
AP_free(msg);
}
}    
AP_check_sends(AP_WAITALL); 
MPI_Barrier(MPI_COMM_WORLD);
delete [] sendcounts;
#ifdef DEBUG
return norient;
#else    
return 0;
#endif  
}
*/

}

#endif
