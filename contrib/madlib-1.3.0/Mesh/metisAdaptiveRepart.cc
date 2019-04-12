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
#ifdef _HAVE_PARMETIS_

#include "mpi.h" // has to be included before "parmetis.h"
extern "C" {
#include "parmetis.h"
}
#include "autopack.h"

#include "assert.h"
#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "ParallelUtils.h"

namespace MAd {

  // -------------------------------------------------------------------
  void metisAdaptiveRepart(pMesh mesh,pMeshDataId tagElt) {

    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);   

    int dim = (mesh->tets.empty()) ? 2 : 3;
    int NN =  (dim == 2) ? mesh->nbTriangles : mesh->nbTets;

    // ------------------------------------------------
    // Creation of 'tab' which contains the number of elements for each proc

    int *tab = new int[nproc];
    int sendnum = NN; 
    if(myrank) {
      MPI_Send(&sendnum,1,MPI_INT,0,myrank,MPI_COMM_WORLD);   
    } else {
      MPI_Status status;
      tab[0] = sendnum;
      for(int i=1 ; i<nproc ; i++) {
        MPI_Recv(&tab[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
      }
    }
    MPI_Bcast(&tab[0],nproc,MPI_INT,0,MPI_COMM_WORLD);
  
    // ------------------------------------------------
    // Creation of 'vtxdist' which contains the starting element id on each proc

    int *vtxdist = new int[nproc+1];
    int pos = 0;
    for(int i=0 ; i<=nproc ; i++) {
      vtxdist[i] = pos;
      if(i<nproc) pos += tab[i];
    }
  
    // ------------------------------------------------
    // Global numbering of the elements

    int mypos = 0;
    for(int i=0; i<myrank; i++) mypos += tab[i];

    if(dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;
      int i=0;
      while ((pface = FIter_next(fit))) { 
        pface->iD = i + mypos;
        i++;
      }
      FIter_delete(fit);  
    }
    else if(dim==3){
      RIter rit = M_regionIter(mesh);
      pRegion pr;
      int i=0;
      while ((pr = RIter_next(rit))) { 
        pr->iD = i + mypos;
        i++;
      }
      RIter_delete(rit);  
    }

    // ------------------------------------------------
    // Creation of the graph of the mesh

    // --- determine the neighborhood of elements ---

    pMeshDataId tagAdj = MD_newMeshDataId("AdjGlobal"); 
    if(dim==2)      E_facesAttachId  (mesh,tagAdj);
    else if(dim==3) F_regionsAttachId(mesh,tagAdj);

    // --- count the adjencies of elements ---

    int *xadj = new int[NN + 1]; // xadj[i] = xadj[i-1] + nb adjencies of element 'i'
    xadj[0] = 0;

    int totCount = 0;

    if(dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;
      int i=0;
      while ((pface = FIter_next(fit))) { 
        int nbAdj = 0;
        assert(!pface->deleted);
        pface->iD = i + mypos;
        for(int k=0 ; k<3 ; k++) {
          pEdge pe = F_edge(pface,k);
          void* tmpptr; 
          int is = EN_getDataPtr((pEntity) pe , tagAdj,&tmpptr);
          assert(is);
          std::vector<int> *recup = (std::vector<int> *) tmpptr;
          nbAdj+=(*recup).size();
        }      
        nbAdj -= 3;
        totCount += nbAdj;  
        xadj[i + 1] = xadj[i] + nbAdj;  
        i++;
      }
      FIter_delete(fit);
    }
    else if(dim==3){
      RIter rit = M_regionIter(mesh);
      pRegion pr;
      int i=0;
      while ((pr = RIter_next(rit))) { 
        int nbAdj = 0;
        for(int k=0 ; k<4 ; k++) {
          pFace pface = R_face(pr,k);
          void* tmpptr; 
          int is = EN_getDataPtr((pEntity) pface , tagAdj,&tmpptr);
          assert(is);
          std::vector<int> *recup = (std::vector<int> *) tmpptr;
          nbAdj+=(*recup).size();
        }      
        nbAdj -= 4;
        totCount += nbAdj;  
        xadj[i + 1] = xadj[i] + nbAdj;  
        i++;
      }
      RIter_delete(rit);
    }

    // --- determine the adjencies ---

    int *adjncy = new int[totCount + 1];

    int count = 0;
    if(dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;
      int i=0;
      while ((pface = FIter_next(fit))) { 
        assert(!pface->deleted);
        for(int k=0 ; k<3 ; k++) {
          pEdge pe = F_edge(pface,k);
          void* tmpptr; 
          int is = EN_getDataPtr((pEntity) pe , tagAdj,&tmpptr);
          assert(is);
          std::vector<int> *recup = (std::vector<int> *) tmpptr;
          for(unsigned int j=0 ; j<(*recup).size() ; j++) {
            int iDrecup = (*recup)[j];
            if(iDrecup != pface->iD)
              adjncy[count++] = iDrecup;
          }            
        }
        i++;
        assert(count==(xadj[i]));
      }
      FIter_delete(fit);
    }
    else {
      RIter rit = M_regionIter(mesh);
      pRegion pr;
      int i=0;
      while ((pr = RIter_next(rit))) { 
        assert(!pr->deleted);
        for(int k=0 ; k<4 ; k++) {
          pFace pface = R_face(pr,k);
          void* tmpptr; 
          int is = EN_getDataPtr((pEntity) pface , tagAdj,&tmpptr);
          assert(is);
          std::vector<int> *recup = (std::vector<int> *) tmpptr;
          for(unsigned int j=0 ; j<(*recup).size() ; j++) {
            int iDrecup = (*recup)[j];
            if(iDrecup != pr->iD)
              adjncy[count++] = iDrecup;
          }            
        }
        i++;
        assert(count==(xadj[i]));
      }
      RIter_delete(rit);
    }  
    
    // ------------------------------------------------
    // Call ParMetis to adapt the graph

    int wgtflag = 0;
    int numflag = 0;
    int ncon    = 1;
    float itr   = 1000.0;
    int options[1];  options[0] = 0;
    int edgecut = 0;
    float *tpwgts, ubvec[1];
    tpwgts  = (float *)   malloc((nproc * ncon) * sizeof(float));
    for (int i=0; i<nproc*ncon; i++)  tpwgts[i] = 1.0/(float)(nproc);
    ubvec[0]  = 1.05;
    MPI_Comm comm = MPI_COMM_WORLD;
    int *partitionVector = new int[NN]; // ParMetis output: will contain the destinations of the elements
    ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, NULL, NULL, NULL,&wgtflag,
                               &numflag,&ncon,&nproc,tpwgts, ubvec,&itr, 
                               options,&edgecut, partitionVector,&comm);
  
    /*    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, NULL, NULL,&wgtflag,
          &numflag,&ncon,&nproc,tpwgts, ubvec, 
          options,&edgecut, partitionVector,&comm);
    */

    // ------------------------------------------------
    // Tag elements to be migrated

    if(dim==2) {
      FIter fit = M_faceIter(mesh);
      pFace pface;
      while ((pface = FIter_next(fit))) {
        int num  = pface->iD - mypos;
        int dest = partitionVector[num];
        if(dest == myrank) continue;
        EN_attachDataInt((pEntity) pface ,tagElt, dest + 1);
      }
      FIter_delete(fit);
    }
    else if(dim==3) {
      RIter rit = M_regionIter(mesh);
      pRegion pr;
      while ((pr = RIter_next(rit))) {
        int num  = pr->iD - mypos;
        int dest = partitionVector[num];
        if(dest == myrank) continue;
        EN_attachDataInt((pEntity) pr ,tagElt, dest + 1);
      }
      RIter_delete(rit);
    }

    // ------------------------------------------------
    // Do some cleaning

    delete [] tab;
    delete [] vtxdist;
    delete [] xadj;
    delete [] partitionVector;

    if(dim==2) {
      EIter eit = M_edgeIter(mesh);
      pEdge pe;
      while ((pe = EIter_next(eit))) { 
        void* tmpptr; 
        int is = EN_getDataPtr((pEntity) pe , tagAdj,&tmpptr);
        assert(is);
        EN_deleteData((pEntity) pe,tagAdj);        
      }    
      EIter_delete(eit);  
    } else if(dim==3) {
      FIter fit = M_faceIter(mesh);
      pFace pface;
      while ((pface = FIter_next(fit))) { 
        void* tmpptr; 
        int is = EN_getDataPtr((pEntity) pface , tagAdj,&tmpptr);
        assert(is);
        EN_deleteData((pEntity) pface,tagAdj);        
      }    
      FIter_delete(fit);
    }
    MD_deleteMeshDataId(tagAdj);
  }

  // -------------------------------------------------------------------
}

#endif
#endif

