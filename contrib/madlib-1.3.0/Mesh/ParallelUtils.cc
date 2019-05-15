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
#include "ParallelUtils.h"
#include "assert.h"

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#endif

#ifdef _HAVE_METIS_
extern "C"
{
#include "metis.h"
}
#endif

namespace MAd {

  // -------------------------------------------------------------------
#ifdef _HAVE_METIS_
  void PartitionMesh(pMesh mesh, int nbPart, const char * filename)
  {
    printf("Partitioning the mesh in %d parts\n",nbPart);

    //NN = number of nodes of the graph

    int dim = (mesh->tets.empty()) ? 2 : 3;
    int NN =  (dim == 2) ? mesh->nbTriangles : mesh->nbTets;
    int * partitionTable = new int[NN];
    int * xadj = new int[NN + 2];

    int totCount = 0;

    if(nbPart>1) {
      MDB_ListF::iterator it2 = mesh->triangles.begin();
      MDB_ListT::iterator it3 = mesh->tets.begin();

      xadj[0] = 0;
      for(int i = 0; i < NN; i++) {
        int nbAdj = 0;
        if(dim == 2) {
          MDB_Triangle *t = *it2;
          ++it2;
          t->iD = i;
          nbAdj = (t->e1->numfaces() + t->e2->numfaces() + t->e3->numfaces() - 3);
          totCount += nbAdj;
        }
        else if(dim == 3) {
          MDB_Tet *t = *it3;
          ++it3;
          t->iD = i;
          nbAdj = (t->f1->getNbRegions() + t->f2->getNbRegions() + 
                   t->f3->getNbRegions() + t->f4->getNbRegions() - 4);
          totCount += nbAdj;
        }
        xadj[i + 1] = xadj[i] + nbAdj;
      }

      it2 = mesh->triangles.begin();
      it3 = mesh->tets.begin();

      int *adjncy = new int[totCount + 1];

      int count = 0;

      for(int i = 0; i < NN; i++) {
        if(dim == 2) {
          MDB_Triangle *t = *it2;
          for(int j = 0; j < t->e1->numfaces(); j++) {
            MDB_Triangle *f = t->e1->faces(j);
            if(f != t)
              adjncy[count++] = f->iD;
          }
          for(int j = 0; j < t->e2->numfaces(); j++) {
            MDB_Triangle *f = t->e2->faces(j);
            if(f != t)
              adjncy[count++] = f->iD;
          }
          for(int j = 0; j < t->e3->numfaces(); j++) {
            MDB_Triangle *f = t->e3->faces(j);
            if(f != t)
              adjncy[count++] = f->iD;
          }
          ++it2;
        }
        else if(dim == 3) {
          MDB_Tet *t = *it3;
          MDB_Tet *o = dynamic_cast<MDB_Tet*>(t->f1->opposite_region((MDB_Region*)t));
          if(o)
            adjncy[count++] = o->iD;
          o = dynamic_cast<MDB_Tet*>(t->f2->opposite_region((MDB_Region*)t));
          if(o)
            adjncy[count++] = o->iD;
          o = dynamic_cast<MDB_Tet*>(t->f3->opposite_region((MDB_Region*)t));
          if(o)
            adjncy[count++] = o->iD;
          o = dynamic_cast<MDB_Tet*>(t->f4->opposite_region((MDB_Region*)t));
          if(o)
            adjncy[count++] = o->iD;
          ++it3;
        }
      }

      int wgtflag = 0;
      int numflag = 0;
      int options[4];
      options[0] = 0;
      int edgecut;
      METIS_PartGraphKway(&NN, xadj, adjncy, 0, 0, &wgtflag,
                          &numflag, &nbPart, options, &edgecut, partitionTable);
      delete[]adjncy;

    } else {
      for(int i = 0; i < NN; i++) {
        partitionTable[i] = 0;     
      }
    }

    M_writeMsh(mesh,filename,2,partitionTable);

    delete[]xadj;
  };
#endif

#ifdef PARALLEL

  // -------------------------------------------------------------------
  struct edge_comm
  {
    pVertex p1,p2;                //pointeur in dest
    int     sendID;       
  };
  struct face_comm
  {
    pVertex p1,p2,p3;                //pointeur in dest
    int     sendID;       
  };

  // -------------------------------------------------------------------
  void E_facesAttachId(pMesh mesh,pMeshDataId tagEdge){
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
  
    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
  
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    EIter eit = M_edgeIter(mesh);
    pEdge pe;  
    while ((pe = EIter_next(eit))) {
      pVertex p1 = E_vertex(pe,0);
      pVertex p2 = E_vertex(pe,1);
      assert(p1); assert(p2);
      assert(E_exist(p1,p2)==pe);
      void *temp_ptr1,*temp_ptr2; 
      int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
      int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
      if(!(isInterface1 && isInterface2))continue;

      const std::vector<std::pair<int , pVertex> > *recup1 = (const std::vector<std::pair<int , pVertex> > *) temp_ptr1;
      const std::vector<std::pair<int , pVertex> > *recup2 = (const std::vector<std::pair<int , pVertex> > *) temp_ptr2;
      int size1 = (*recup1).size();
      int size2 = (*recup2).size();
      assert(size1);assert(size2);
      int *tab=new int[size1];
      for(int i=0 ; i< size1 ; i++) tab[i] = (*recup1)[i].first;
      for(int j=0 ; j<size2 ; j++) {    
        int iProc = (*recup2)[j].first;
        int i;
        for(i=0 ; i<size1 ; i++) {
          if(iProc == tab[i]) break;
        }
        if(i < size1) {        
          pVertex remote1 = (*recup1)[i].second;
          pVertex remote2 = (*recup2)[j].second;
          for(int k = 0; k < pe->numfaces(); k++) {
            pFace pface = pe->faces(k);
            assert(pface);
            assert(iProc != myrank);
            void *buf = AP_alloc(iProc,444,sizeof(edge_comm));
            edge_comm *castbuf = (edge_comm *) buf;
            castbuf->p1	     = remote1;
            castbuf->p2	     = remote2;
            castbuf->sendID    = pface->iD;
            AP_send(buf);
            sendcounts[iProc]++;  
          }
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
        int numrecv = castbuf->sendID;
        int NN = mesh->nbTriangles;
        if(from > myrank) {
          if(numrecv < NN) 
            printf("from %d iD %d NN %d \n",from,numrecv,NN);
          assert(numrecv >= NN);
        }
        assert(from!=myrank);
        pEdge pe = E_exist(p1recv,p2recv);
        if(pe) {
          void* tmpptr; 
          int is = EN_getDataPtr((pEntity) pe , tagEdge,&tmpptr);
          if(is){
            std::vector<int> *recup = (std::vector<int> *) tmpptr;
            for(unsigned int i=0 ; i<(*recup).size() ; i++) assert((*recup)[i]!=numrecv);	     
            (*recup).push_back(numrecv);
          } else {
            std::vector<int> list;
            list.push_back(numrecv);
            EN_attachDataPtr((pEntity) pe,tagEdge,new std::vector<int>(list));       
          }
        }
        AP_free(msg);
      }
    }    
    AP_check_sends(AP_WAITALL); 
    MPI_Barrier(MPI_COMM_WORLD);

    delete [] sendcounts;
  
    eit = M_edgeIter(mesh);
    while ((pe = EIter_next(eit))) { 
      void* tmpptr; 
      int is = EN_getDataPtr((pEntity) pe , tagEdge,&tmpptr);
      if(is){
        std::vector<int> *recup = (std::vector<int> *) tmpptr;
        for(int k = 0; k <  E_numFaces(pe); k++) {
          pFace pface = E_face(pe,k);
          int iD      = pface->iD;
          unsigned int i;
          for(i=0 ; i<(*recup).size() ; i++) {
            if((*recup)[i]==iD) {
              printf("%d my %d iD %d \n",myrank,(*recup)[i],pe->numfaces());
            }
            assert((*recup)[i]!=iD);
          }  
          (*recup).push_back(iD);
        }
      } else {
        std::vector<int> list;
        assert(pe->numfaces());
        for(int k = 0; k < pe->numfaces(); k++) {
          pFace pface = pe->faces(k);
          list.push_back(pface->iD);
        }	 
        EN_attachDataPtr((pEntity) pe,tagEdge,new std::vector<int>(list));        
      }    
    } 
    EIter_delete(eit);

  }

  // -------------------------------------------------------------------
  void F_regionsAttachId(pMesh mesh,pMeshDataId tagAdj){
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
  
    int *sendcounts = new int[nproc];
    for(int i=0;i<nproc;i++)sendcounts[i]=0;
  
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    FIter fit = M_faceIter(mesh);
    pFace pface;  
    while ((pface = FIter_next(fit))) {
      pVertex p1 = F_vertex(pface,0);
      pVertex p2 = F_vertex(pface,1);
      pVertex p3 = F_vertex(pface,2);
      assert(p1); assert(p2);assert(p3);
      void *temp_ptr1,*temp_ptr2,*temp_ptr3; 
      int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
      int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
      int isInterface3 = EN_getDataPtr((pEntity) p3 , tagData, &temp_ptr3);
      if(!(isInterface1 && isInterface2 && isInterface3))continue;

      const std::vector<std::pair<int , pVertex> > *recup1 = (const std::vector<std::pair<int , pVertex> > *) temp_ptr1;
      const std::vector<std::pair<int , pVertex> > *recup2 = (const std::vector<std::pair<int , pVertex> > *) temp_ptr2;
      const std::vector<std::pair<int , pVertex> > *recup3 = (const std::vector<std::pair<int , pVertex> > *) temp_ptr3;
      int size1 = (*recup1).size();
      int size2 = (*recup2).size();
      int size3 = (*recup3).size();
      assert(size1);assert(size2);assert(size3);
      int *tab=new int[size1];
      int *tab3=new int[size3];
      for(int i=0 ; i< size1 ; i++) tab[i]  = (*recup1)[i].first;
      for(int i=0 ; i< size3 ; i++) tab3[i] = (*recup3)[i].first;
      for(int j=0 ; j<size2 ; j++) {    
        int iProc = (*recup2)[j].first;
        int i,l;
        for(i=0 ; i<size1 ; i++) {
          if(iProc == tab[i]) break;
        }
        for(l=0 ; l<size3 ; l++) {
          if(iProc == tab3[l]) break;
        }
        if(i < size1 && l<size3) {        
          pVertex remote1 = (*recup1)[i].second;
          pVertex remote2 = (*recup2)[j].second;
          pVertex remote3 = (*recup3)[l].second;
          for(int k = 0; k < F_numRegions(pface); k++) {
            pRegion pr = F_region(pface,k);
            assert(pr);
            assert(iProc != myrank);
            void *buf = AP_alloc(iProc,444,sizeof(face_comm));
            face_comm *castbuf = (face_comm *) buf;
            castbuf->p1	     = remote1;
            castbuf->p2	     = remote2;
            castbuf->p3	     = remote3;
            castbuf->sendID  = pr->iD;
            AP_send(buf);
            sendcounts[iProc]++;  
          }
        }
      }
      delete []tab;
      delete []tab3;
    }      
    FIter_delete(fit);  
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
        face_comm * castbuf = (face_comm*) msg;
        pVertex p1recv,p2recv,p3recv;
        p1recv = castbuf -> p1;
        p2recv = castbuf -> p2;
        p3recv = castbuf -> p3;
        int numrecv = castbuf->sendID;
        int NN = mesh->nbTets;
        if(from > myrank) {
          if(numrecv < NN) 
            printf("from %d iD %d NN %d \n",from,numrecv,NN);
          assert(numrecv >= NN);
        }
        assert(from!=myrank);
        // pFace pface = F_exist(2,p1recv,p2recv,p3recv,0);
        pFace pface = F_exist(p1recv,p2recv,p3recv,0);
        if(pface) {
          void* tmpptr; 
          int is = EN_getDataPtr((pEntity) pface , tagAdj,&tmpptr);
          if(is){
            std::vector<int> *recup = (std::vector<int> *) tmpptr;
            for(unsigned int i=0 ; i<(*recup).size() ; i++) assert((*recup)[i]!=numrecv);	     
            (*recup).push_back(numrecv);
          } else {
            std::vector<int> newvec;
            newvec.push_back(numrecv);
            EN_attachDataPtr((pEntity) pface,tagAdj,new std::vector<int>(newvec));       
          }
        }
        AP_free(msg);
      }
    }    
    AP_check_sends(AP_WAITALL); 
    MPI_Barrier(MPI_COMM_WORLD);

    delete [] sendcounts;
  
    fit = M_faceIter(mesh);
    while ((pface = FIter_next(fit))) { 
      void* tmpptr; 
      int is = EN_getDataPtr((pEntity) pface , tagAdj,&tmpptr);
      if(is){
        std::vector<int> *recup = (std::vector<int> *) tmpptr;
        for(int k = 0; k <  F_numRegions(pface); k++) {
          pRegion pr = F_region(pface,k);
          int iD      = pr->iD;
          unsigned int i;
          for(i=0 ; i<(*recup).size() ; i++) {
            if((*recup)[i]==iD) {
              printf("%d my %d iD %d \n",myrank,(*recup)[i],F_numRegions(pface));
            }
            assert((*recup)[i]!=iD);
          }  
          (*recup).push_back(iD);
        }
      } else {
        std::vector<int> list;
        assert(F_numRegions(pface));
        for(int k = 0; k < F_numRegions(pface); k++) {
          pRegion pr = F_region(pface,k);
          list.push_back(pr->iD);
        }	 
        EN_attachDataPtr((pEntity) pface,tagAdj,new std::vector<int>(list));        
      }    
    } 
    FIter_delete(fit);

  }

  // -------------------------------------------------------------------
#endif

}
