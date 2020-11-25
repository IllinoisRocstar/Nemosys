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

#include "Mark.h"
#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"

#include "assert.h"

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#endif

namespace MAd {

#ifdef PARALLEL
  // -------------------------------------------------------------------
  void MarkTets(pMesh mesh,pMeshDataId tagElt)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");

    int * tab = new int[nproc];
    int sendnum = mesh->nbPoints; 
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
    MPI_Barrier(MPI_COMM_WORLD);

    int npGlob = 0;
    for(int i=0 ; i<nproc ; i++) npGlob +=tab[i];
  
    /* 1) marked regions*/
    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      //int dest  = rand()%nproc;

      int dest = myrank;// = rand()%nproc;
      pVertex nod[4];
      for(int i=0 ; i<4 ; i++) {
        nod[i] = R_vertex(pr,i);
        void *temp_ptr; 
        int isInterface = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
        int maxproc = myrank;
        if(isInterface) {
          const std::vector<std::pair<int , pVertex> > *recup = (std::vector<std::pair<int , pVertex> > *) temp_ptr;
          for(unsigned int j=0 ; j<(*recup).size() ; j++) {
            int numproc = (*recup)[j].first;
            if(tab[numproc] < tab[maxproc]) {
              maxproc = numproc;
            }
          }
        }
        if(tab[maxproc] < tab[dest]) dest = maxproc;
      }  
      assert(dest<nproc);
      assert(dest>=0);
      if(dest == myrank) continue;
      EN_attachDataInt((pEntity) pr ,tagElt, dest + 1);
    }  
    RIter_delete(rit);
    delete []tab;
  }


  // -------------------------------------------------------------------
  void MarkTetsSmooth(pMesh mesh,pMeshDataId tagElt)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");

    int *tab=new int[nproc];
    int sendnum = mesh->nbPoints; 
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
    MPI_Barrier(MPI_COMM_WORLD);

    int npGlob = 0;
    for(int i=0 ; i<nproc ; i++) npGlob +=tab[i];
  
    /* 1) marked regions*/
    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      //int dest  = rand()%nproc;

      int dest = myrank;// = rand()%nproc;
      pVertex nod[4];
      for(int i=0 ; i<4 ; i++) {
        nod[i] = R_vertex(pr,i);
        void *temp_ptr; 
        int isInterface = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
        int maxproc = myrank;
        if(isInterface) {
          const std::vector<std::pair<int , pVertex> > *recup = (std::vector<std::pair<int , pVertex> > *) temp_ptr;
          for(unsigned int j=0 ; j<(*recup).size() ; j++) {
            int numproc = (*recup)[j].first;
            if(tab[numproc] < tab[maxproc]) {
              maxproc = numproc;
            }
          }
        }
        if(tab[maxproc] < tab[dest]) dest = maxproc;
      }  
      assert(dest<nproc);
      assert(dest>=0);
      if(dest == myrank) continue;
      EN_attachDataInt((pEntity) pr ,tagElt, dest + 1);
    }  
    RIter_delete(rit);

    //mark neighboor of interface tetras
    int    iter = 0;
    int maxiter = 1;
    int  tetmove = 10;  
    while(iter<maxiter) {
      tetmove = 0;
      rit = M_regionIter(mesh); 
      while ((pr = RIter_next(rit))) {
        int dest;
        int isMove = EN_getDataInt((pEntity) pr ,tagElt, &dest);
        if(!isMove) continue;
        if(iter%2)
          if(dest < 0) continue;
          else 
            if(dest > 0) continue;
        EN_attachDataInt((pEntity) pr ,tagElt, -dest) ;
        for(int i = 0 ; i<4 ; i++) {
          pFace pface     = R_face(pr,i);
          pRegion prother;
          int k =0;
          for(k = 0; k <  F_numRegions(pface); k++) {
            prother = F_region(pface,k);
            if(prother != pr) break;
          }
          if(k ==	F_numRegions(pface)) continue;
          int destother;
          int is = EN_getDataInt((pEntity) prother ,tagElt, &destother);
          if(!is) EN_attachDataInt((pEntity) prother ,tagElt, -dest) ;    
        }
      }  
      RIter_delete(rit);  
      //printf("iter neigh %d : %d\n",iter,trmove);
      iter++;
    }
    //smooth new interfaces
    iter = 0;
    maxiter = 10;
    tetmove = 10;
    while(tetmove && iter<maxiter) {
      tetmove = 0;
      rit = M_regionIter(mesh); 
      while ((pr = RIter_next(rit))) {
        int dest;
        int isMove = EN_getDataInt((pEntity) pr ,tagElt, &dest);
        if(!isMove) continue;
        if(dest < 0) EN_attachDataInt((pEntity) pr ,tagElt, abs(dest));
        dest = abs(dest);
        if((dest-1)== myrank) continue;
        //how many neighboors in dest-1?
        int nbN = 0;
        int destreal = dest;
        int nb = 0;
        for(int i = 0 ; i < 4 ; i++) {
          pFace pface     = R_face(pr,i);
          pRegion prother;
          int k =0;
          for(k = 0; k <  F_numRegions(pface); k++) {
            prother = F_region(pface,k);
            if(prother != pr) break;
          }
          if(k ==	F_numRegions(pface)) continue;
          int destother;
          int is = EN_getDataInt((pEntity) prother ,tagElt, &destother);
          if(!is) continue;
          if(destother == dest) nb++;
          else destreal = destother;    
        }
        int nbGood = 4 - nbN + nb;
        if(nbGood==1 || nbGood==2) {
          tetmove++;
          //printf("smoothing new interfaces : %d\n",pface->iD);
          if(destreal == dest) {
            //printf("nbGood : %d %d\n",nbN,nb);
            EN_deleteData((pEntity) pr ,tagElt);
          } else {
            EN_attachDataInt((pEntity) pr ,tagElt, abs(destreal));
          }
        }
      }  
      RIter_delete(rit);
      //printf("iter smoothing %d : %d\n",iter,tetmove);
      iter++;
    }
  
    //attach positive destination
    rit = M_regionIter(mesh);
    while ((pr = RIter_next(rit))) {
      int dest;
      int isMove = EN_getDataInt((pEntity) pr ,tagElt, &dest);
      if(isMove) {
        if(dest < 0) EN_attachDataInt((pEntity) pr ,tagElt, abs(dest));
        //assert(dest > 0);
      }
    }  
    RIter_delete(rit);
    delete []tab;
  }

  // -------------------------------------------------------------------
  void MarkTetsRandom(pMesh mesh, pMeshDataId tagElt) {

    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      int dest = rand() % nproc;
      if(dest == myrank) continue;
      EN_attachDataInt((pEntity) pr, tagElt, dest + 1);
    }  
    RIter_delete(rit);
  }

  // -------------------------------------------------------------------
  int MarkTetsManifold(pMesh mesh,pMeshDataId tagElt) {

    int nmanifold = 0;

    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    pMeshDataId tagData   = MD_lookupMeshDataId("RemotePoint");
    EIter eit = M_edgeIter(mesh);
    pEdge ped;
    while ((ped = EIter_next(eit))) {
      pVertex p1 = E_vertex(ped,0);
      pVertex p2 = E_vertex(ped,1);
      void *temp_ptr1,*temp_ptr2; 
      int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
      int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
      if(!(isInterface1 && isInterface2)) continue;
      int nFace = 0;
      int dest = myrank;// = rand()%nproc;
      int maxproc = myrank;
      for(int i=0 ; i<E_numFaces(ped) ; i++) {
        pFace pface = E_face(ped,i);
        for(int j=0 ; j<F_numRegions(pface) ; j++) {
          pRegion pr = F_region(pface,j);
          pVertex nod[4];
          for(int i=0 ; i<4 ; i++) {
            nod[i] = R_vertex(pr,i);
            void *temp_ptr; 
            int isInterface = EN_getDataPtr((pEntity) nod[i] , tagData, &temp_ptr);
            int maxproc = myrank;
            if(isInterface) {
              const std::vector<std::pair<int , pVertex> > *recup = (std::vector<std::pair<int , pVertex> > *) temp_ptr;
              for(unsigned int j=0 ; j<(*recup).size() ; j++) {
                int numproc = (*recup)[j].first;
                if(numproc < maxproc) {
                  maxproc = numproc;
                }
              }
            }
          }
          if(maxproc < dest) dest = maxproc;  	
        }
        if(F_numRegions(pface)==2) continue;
        int j;
        for(j=0 ; j<3 ; j++) {
          pVertex pp = F_vertex(pface,j);
          void *temp_ptr; 
          int isInterface = EN_getDataPtr((pEntity) pp , tagData, &temp_ptr);
          if(!isInterface) break;
        }
        if(j==3) nFace++;
      }
      if(nFace<=2) continue;
      nmanifold++;
      //printf("nFace %d non manifold edge\n",nFace);
      if(dest == myrank) continue;
      for(int i=0 ; i<E_numFaces(ped) ; i++) {
        pFace pface = E_face(ped,i);
        for(int j=0 ; j<F_numRegions(pface) ; j++) {
          pRegion pr = F_region(pface,j);
          EN_attachDataInt((pEntity) pr ,tagElt, dest + 1);
        }
      }
    }
    EIter_delete(eit);
  
    return(nmanifold);  
  }

  // -------------------------------------------------------------------
  void MarkTriangles(pMesh mesh, pMeshDataId tagElt)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");

    // --- Every proc send its number of nodes -> collected in tab ---
    int * tab = new int[nproc];
    int sendnum = M_numVertices(mesh); 
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
    MPI_Barrier(MPI_COMM_WORLD);

    int npGlob = 0;
    for(int i=0 ; i<nproc ; i++) npGlob +=tab[i];
  
    // --- Mark the triangles which have at least one vertex on an interface 
    //     with a distant mesh that has less nodes than the current one ---
    FIter fit = M_faceIter(mesh);
    pFace pface;
    while ( (pface = FIter_next(fit)) ) {
      pVertex nod[3];
      pface->getNodes(nod);
      int dest = myrank;
      for(int i=0; i<3; i++) {
        void *temp_ptr;
        int isInterface = EN_getDataPtr((pEntity) nod[i], tagData, &temp_ptr);
        int maxproc = myrank;
        if(isInterface) {
          const std::vector<std::pair<int, pVertex> > *recup = (const std::vector<std::pair<int, pVertex> > *) temp_ptr;
          for(int j=0; j<(int)((*recup).size()); j++) {
            int numproc = (*recup)[j].first;
            if(tab[numproc] < tab[maxproc]) maxproc = numproc;
          }
        }
        if(tab[maxproc] < tab[dest]) dest = maxproc;
      }  
      assert(dest<nproc);
      assert(dest>=0);
      if(dest == myrank) continue;
      EN_attachDataInt((pEntity) pface, tagElt, dest + 1);
    }  
    FIter_delete(fit);
    delete []tab;
  }


  
  // -------------------------------------------------------------------
  void MarkTrianglesSmooth(pMesh mesh,pMeshDataId tagElt)
  {
    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");

    // --- Every proc send its number of nodes -> collected in tab ---
    int * tab = new int[nproc];
    int sendnum = M_numVertices(mesh); 
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
    MPI_Barrier(MPI_COMM_WORLD);

    int npGlob = 0;
    for(int i=0 ; i<nproc ; i++) npGlob += tab[i];
  
    // --- Mark the triangles which have at least one vertex on an interface 
    //     with a distant mesh that has less nodes than the current one ---
    FIter fit = M_faceIter(mesh);
    pFace pface;
    while ( (pface = FIter_next(fit)) ) {
      pVertex nod[3];
      pface->getNodes(nod);
      int dest = myrank;
      for(int i=0; i<3; i++) {
        void *temp_ptr;
        int isInterface = EN_getDataPtr((pEntity) nod[i], tagData, &temp_ptr);
        int maxproc = myrank;
        if(isInterface) {
          const std::vector<std::pair<int, pVertex> > *recup = (const std::vector<std::pair<int, pVertex> > *) temp_ptr;
          for(int j=0; j<(int)((*recup).size()); j++) {
            int numproc = (*recup)[j].first;
            if(tab[numproc] < tab[maxproc]) maxproc = numproc;
          }
        }
        if(tab[maxproc] < tab[dest]) dest = maxproc;
      }  
      assert(dest<nproc);
      assert(dest>=0);
      if(dest == myrank) continue;
      EN_attachDataInt((pEntity) pface, tagElt, dest + 1);
    }  
    FIter_delete(fit);
  
    // --- Mark the neighbour triangles of the marked triangles (do it n times) ---
    // GCREMARK: why not apply the same principle as for the first marking -> mark around nodes ?
    int    iter = 0;
    int maxiter = 1;
    int  trmove = 10;  
    while(iter<maxiter) {
      trmove = 0;
      fit = M_faceIter(mesh); 
      while ((pface = FIter_next(fit))) {
        int dest;
        int isMove = EN_getDataInt((pEntity) pface ,tagElt, &dest);
        if(!isMove) continue;
        if(iter%2)
          if(dest < 0) continue;
          else 
            if(dest > 0) continue;
        EN_attachDataInt((pEntity) pface ,tagElt, -dest) ;
        for(int i = 0 ; i<3 ; i++) {
          pEdge ped     = F_edge(pface,i);
          pFace pfother;
          int k =0;
          for(k = 0; k <  E_numFaces(ped); k++) {
            pfother = E_face(ped,k);
            if(pfother != pface) break;
          }
          if(k ==	E_numFaces(ped)) continue;
          int destother;
          int is = EN_getDataInt((pEntity) pfother ,tagElt, &destother);
          if(!is) EN_attachDataInt((pEntity) pfother ,tagElt, -dest) ;    
        }
      }  
      FIter_delete(fit);  
      //printf("iter neigh %d : %d\n",iter,trmove);
      iter++;
    }

    // --- Smooth new interfaces ---
    iter = 0;
    maxiter = 10;
    trmove = 10;
    while(trmove && iter<maxiter) {
      trmove = 0;
      fit = M_faceIter(mesh); 
      while ((pface = FIter_next(fit))) {
        int dest;
        int isMove = EN_getDataInt((pEntity) pface ,tagElt, &dest);
        if(!isMove) continue;
        if(dest < 0) EN_attachDataInt((pEntity) pface ,tagElt, abs(dest));
        dest = abs(dest);
        if((dest-1)== myrank) continue;
        //how many neighbours in dest-1?
        int nbN = 0;
        int destreal = dest;
        int nb = 0;
        for(int i = 0 ; i < 3 ; i++) {
          pEdge ped     = F_edge(pface,i);
          pFace pfother;
          int k =0;
          for(k = 0; k <  E_numFaces(ped); k++) {
            pfother = E_face(ped,k);
            if(pfother != pface) break;
          }
          if(k ==	E_numFaces(ped)) continue;
          int destother;
          int is = EN_getDataInt((pEntity) pfother ,tagElt, &destother);
          if(!is) continue;
          if(destother == dest) nb++;
          else destreal = destother;    
        }
        int nbGood = 3 - nbN + nb;
        if(nbGood==1) {
          trmove++;
          //printf("smoothing new interfaces : %d\n",pface->iD);
          if(destreal == dest) {
            //printf("nbGood : %d %d\n",nbN,nb);
            EN_deleteData((pEntity) pface ,tagElt);
          } else {
            EN_attachDataInt((pEntity) pface ,tagElt, abs(destreal));
          }
        }
      }  
      FIter_delete(fit);
      printf("iter smoothing %d : %d\n",iter,trmove);
      iter++;
    }

    // --- Attach positive destination ---
    fit = M_faceIter(mesh);
    while ((pface = FIter_next(fit))) {
      int dest;
      int isMove = EN_getDataInt((pEntity) pface ,tagElt, &dest);
      if(isMove) {
        if(dest < 0) EN_attachDataInt((pEntity) pface ,tagElt, abs(dest));
        //assert(dest > 0);
      }
    }  
    FIter_delete(fit);
    delete []tab;
  }

  // -------------------------------------------------------------------
  void MarkTrianglesRandom(pMesh mesh, pMeshDataId tagElt) {

    int nproc,myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
    FIter fit = M_faceIter(mesh);
    pFace pf;
    while ((pf = FIter_next(fit))) {
      int dest = rand() % nproc;
      if(dest == myrank) continue;
      EN_attachDataInt((pEntity) pf, tagElt, dest + 1);
    }
    FIter_delete(fit);
  }
#endif

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  /*loop on an edges list for 3-periodic cases*/   
  bool compare(pEdge p1, pEdge p2) {
    pVertex p10 = E_vertex(p1,0);
    pVertex p11 = E_vertex(p1,1);
    pVertex p20 = E_vertex(p2,0);
    pVertex p21 = E_vertex(p2,1);

    double len1 = (p10->X-p11->X)*(p10->X-p11->X) + (p10->Y-p11->Y)*(p10->Y-p11->Y) +(p10->Z-p11->Z)*(p10->Z-p11->Z);
    double len2 = (p20->X-p21->X)*(p20->X-p21->X) + (p20->Y-p21->Y)*(p20->Y-p21->Y) +(p20->Z-p21->Z)*(p20->Z-p21->Z);
    //	printf("    edge %d %d : %e\n",EN_id((pEntity) p10),EN_id((pEntity) p11),len1) ;
    //	printf(" et edge %d %d : %e\n",EN_id((pEntity) p20),EN_id((pEntity) p21),len2) ;
    int num11 = EN_id((pEntity) p10);
    int num12 = EN_id((pEntity) p11);
    int num21 = EN_id((pEntity) p20);
    int num22 = EN_id((pEntity) p21);  
    //	printf("egal %d %d %d\n",len1==len2,num11==num21,num12<num22)  ;
	
    if(len1==len2) {
      if(num11==num21) {
        return(num12>num22);
      } else {
        return(num11>num21);  
      }
    }
    return ((len1 > len2));
  }  

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void VertexTagMoveList(pMesh mesh, std::vector<std::vector<int> >& transforef,
                         pMeshDataId tagMove,pMeshDataId tagTransfo) {
    pMeshDataId tagPeriodic     = MD_lookupMeshDataId("PeriodicPoint");
  
    /*building list*/ 
    std::list<pEdge> listedge;
    EIter eit = M_edgeIter(mesh);
    pEdge ped;  
    while ((ped = EIter_next(eit))) { 
      //	pVertex p0 = E_vertex(ped,0);
      //    pVertex p1 = E_vertex(ped,1);
      //printf("adding %d %d\n",EN_id((pEntity) p0),EN_id((pEntity) p1));
      listedge.push_back(ped);
    } 
    EIter_delete(eit);  
    //GCRMK: anisotropic case!!
    listedge.sort(compare);   

    /*favorite direction ?*/
    int direct = 0;
    int refsize = transforef.size();
    for(unsigned int kt = 0 ; kt < transforef[0].size(); kt++) {
      direct += transforef[0][kt]*transforef[0][kt];
    }
    //printf("********** favorite direction ??? %d\n",direct);


    for(std::list<pEdge>::iterator it = listedge.begin() ; it!=listedge.end() ; it++) {  
      pEdge ped = (*it);
      /*summits of edge ped*/
      pVertex p0 = E_vertex(ped,0);
      pVertex p1 = E_vertex(ped,1);
      if(it==listedge.begin()) printf("the longest is %d %d\n",EN_id(p0),EN_id(p1));
      pVertex pFix0/*,pFix1*/;
      std::vector<int> vectnod;
      int find = 0; 	
      /*sommets periodics ?*/
      void *temp_ptr0; 
      int isP0 = EN_getDataPtr((pEntity) p0 , tagPeriodic, &temp_ptr0);
      if(!isP0) continue;
      void *temp_ptr1; 
      int isP1 = EN_getDataPtr((pEntity) p1 , tagPeriodic, &temp_ptr1);
      if(!isP1) continue; 
      /*sommets deja tagges ?*/
      int move0;
      int isM0 = EN_getDataInt((pEntity) p0 ,tagMove, &move0);
      if(!(!isM0 || (move0 == -2))) continue;	  
      int move1;
      int isM1 = EN_getDataInt((pEntity) p1 ,tagMove, &move1);
      if(!(!isM1 || (move1 == -2))) continue;	  
    
      //printf("treat edge %d %d\n",EN_id((pEntity) p0),EN_id((pEntity) p1));    
      std::vector<std::pair<std::vector<int> , pVertex> > *recup0 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr0;
      std::vector<std::pair<std::vector<int> , pVertex> > *recup1 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;
    
      /*a-t-on deja une direction privilegiee ?*/
      if(1/*!direct*/) {
        unsigned int j=0;
        for(j=0 ; j<(*recup0).size() ; j++) {
          std::vector<int> transfo0 = (*recup0)[j].first;
          pVertex pImg0 = (*recup0)[j].second;
          //if(pImg0==p1) continue; 
          unsigned int j1;
          for(j1=0 ; j1<(*recup1).size() ; j1++) {
            std::vector<int> transfo1 = (*recup1)[j1].first;
            pVertex pImg1 = (*recup1)[j1].second;
            //if(pImg1==p0) continue;  
    
            /*transfo0 == transfo1 ?*/
            assert(transfo1.size()==transfo0.size());
            unsigned int kt = 0;
            for(kt = 0 ; kt < transfo0.size(); kt++) {
              if(transfo0[kt] != transfo1[kt]) break;
            }
            if(kt!=transfo0.size()) continue; 
            /* pImg1-pImg0 existe ?*/
            if ( !E_exist(pImg0,pImg1) && !E_exist(pImg1,pImg0) ) continue;
            /* a-t-on deja trouve une arete img ?*/
            //printf("img edge %d %d\n",EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));    
            if(!find) {
              //if(it==listedge.begin()) printf("c'est bon on bouge\n");
              /*si le point img bouge, on bouge pas!*/
              int moveI0;
              int isMI0 = EN_getDataInt((pEntity) pImg0 ,tagMove, &moveI0);
              int moveI1;
              int isMI1 = EN_getDataInt((pEntity) pImg1 ,tagMove, &moveI1);
              if (isMI0 && (moveI0==1) ) continue;   	  
              if (isMI1 && (moveI1==1) ) continue;  
              if(!direct) { 	  
                find = 1;
                int move = 1;  
                EN_attachDataInt((pEntity) p0 ,tagMove, move);
                EN_attachDataInt((pEntity) p1 ,tagMove, move);
                move = -1;  
                EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
                EN_attachDataInt((pEntity) pImg1 ,tagMove, move);
    
                /*on garde les deux points img*/
                pFix0 = pImg0;
                //pFix1 = pImg1;
                /*rajout de la transfo avec laquelle bouger les points*/
                for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
                  vectnod.push_back(transfo0[kt]);
                  transforef[0][kt] = transfo0[kt];         
                } 
                EN_attachDataPtr((pEntity) p0 , tagTransfo,new std::vector<int>(vectnod)); 
                EN_attachDataPtr((pEntity) p1 , tagTransfo,new std::vector<int>(vectnod)); 	
                for(unsigned int kt = 0 ; kt < transforef[0].size(); kt++) {
                  direct += transforef[0][kt]*transforef[0][kt];
                } 
                //printf("on a une direct %d : %d %d\n",direct,transforef.size(),refsize);
              } else { /*on a deja une direction privilegiee*/
                double pds = 0;
                for(int st = 0 ; st < refsize; st++) {
                  for(unsigned int kt = 0 ; kt < transforef[st].size(); kt++) {
                    pds += transforef[st][kt]*transfo0[kt];
                  }
                  if(pds!=0) break; 
                }
                if(pds>=0) {
                  //printf("direction ok : %d %d %d (%d)\n",transfo0[0],transfo0[1],transfo0[2],transforef.size()); 
                  if(pds==0) {/*on rajoute une direction privilegiee*/
                    transforef.push_back(transfo0);
                    refsize++;       
                    //printf("on rajoute la direction au ref : %d \n",refsize);
                  }   
                  find = 1;
                  int move = 1;  
                  EN_attachDataInt((pEntity) p0 ,tagMove, move);
                  EN_attachDataInt((pEntity) p1 ,tagMove, move);
                  move = -1;  
                  EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
                  EN_attachDataInt((pEntity) pImg1 ,tagMove, move);
    
                  /*on garde les deux points img*/
                  pFix0 = pImg0;
                  //pFix1 = pImg1;
                  /*rajout de la transfo avec laquelle bouger les points*/
                  for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
                    vectnod.push_back(transfo0[kt]);         
                  }
                  EN_attachDataPtr((pEntity) p0 , tagTransfo,new std::vector<int>(vectnod)); 
                  EN_attachDataPtr((pEntity) p1 , tagTransfo,new std::vector<int>(vectnod)); 	
                } else {/*la direction est opposee a celle souhaitee --> on fait l'inverse*/
                  //printf("bad direction : %d %d %d (%d %d)\n",transfo0[0],transfo0[1],transfo0[2],transforef.size(),refsize);    
                  /*****INUTILE ????????******/
                  int moveI0;
                  int isMI0 = EN_getDataInt((pEntity) p0 ,tagMove, &moveI0);
                  int moveI1;
                  int isMI1 = EN_getDataInt((pEntity) p1 ,tagMove, &moveI1);
                  if (isMI0 && (moveI0==1) ) break;  	  
                  if (isMI1 && (moveI1==1) ) break;     	  
                  if(!(!isMI0 || (moveI0 == -2))) break;     
                  if(!(!isMI0 || (moveI0 == -2))) break;  
                  /*****FIN INUTILE ????????******/
                  find = 1;
                  int move = -1;  
                  EN_attachDataInt((pEntity) p0 ,tagMove, move);
                  EN_attachDataInt((pEntity) p1 ,tagMove, move);
                  move = 1;  
                  EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
                  EN_attachDataInt((pEntity) pImg1 ,tagMove, move);

                  /*on garde les deux points img*/
                  pFix0 = p0;
                  //pFix1 = p1;
                  /*rajout de la transfo avec laquelle bouger les points*/
                  for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
                    vectnod.push_back((-1)*transfo0[kt]);         
                  }
                  EN_attachDataPtr((pEntity) pImg0 , tagTransfo,new std::vector<int>(vectnod)); 
                  EN_attachDataPtr((pEntity) pImg1 , tagTransfo,new std::vector<int>(vectnod)); 	
                }  
              }  /*end if direct*/        
              break;
            } else {
              /*si le point bouge deja, on veut pas chger sa destination!!!*/
              int moveI0;
              int isMI0 = EN_getDataInt((pEntity) pImg0 ,tagMove, &moveI0);
              int moveI1;
              int isMI1 = EN_getDataInt((pEntity) pImg1 ,tagMove, &moveI1);
              if (isMI0 && (moveI0==1) ) continue;   	  
              if (isMI1 && (moveI1==1) ) continue;   	  
          
              int move = 1; 
              EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
              EN_attachDataInt((pEntity) pImg1 ,tagMove, move); 
          
              /*on sait que pImg0pImg1 doit aller sur pFix0pFix1 il faut trouver la transfo*/
              /*si PFix0==p0   transfo = -transfo0*/ 
              /*sinon          transfo = -transfo0 + vecnod  */
              std::vector<int> vecttrans;
              if(pFix0==p0) {
                for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
                  vecttrans.push_back(-transfo0[kt]);
                }			
              } else {
                for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
                  vecttrans.push_back(-transfo0[kt]);
                }
                for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
                  vecttrans[kt] += vectnod[kt];
                }			
              }     
              EN_attachDataPtr((pEntity) pImg0 , tagTransfo,new std::vector<int>(vecttrans)); 
              EN_attachDataPtr((pEntity) pImg1 , tagTransfo,new std::vector<int>(vecttrans)); 		      
              break;			
            }	    
          } /*end j1*/
          if(j1==(*recup1).size()){
            //printf("arete not found\n");
            /* int move = -2;  
               EN_attachDataInt((pEntity) pImg0 ,tagMove, move);  */  	  
          }
        } /*end j*/	  
  	  //if(!find) {printf("on n'a pas trouve d'arete img donc on tag pas\n");       exit(0);}
      } else {/*if direct*/
	
      }    
    }       
    return;
  }


  // -------------------------------------------------------------------
  /*boucle sur les aretes pour traiter les cas 3-periodic*/
  void VertexTagMove(pMesh mesh,pMeshDataId tagMove,pMeshDataId tagTransfo) {
    pMeshDataId tagPeriodic     = MD_lookupMeshDataId("PeriodicPoint");

    EIter eit = M_edgeIter(mesh);
    pEdge ped,pedla;  
    double lenla = 0;
    while ((ped = EIter_next(eit))) { 
      /*les deux sommets de l'arete ped*/
      pVertex p0 = E_vertex(ped,0);
      pVertex p1 = E_vertex(ped,1);

      /*sommets periodics ?*/
      void *temp_ptr0; 
      int isP0 = EN_getDataPtr((pEntity) p0 , tagPeriodic, &temp_ptr0);
      if(!isP0) continue;
      void *temp_ptr1; 
      int isP1 = EN_getDataPtr((pEntity) p1 , tagPeriodic, &temp_ptr1);
      if(!isP1) continue; 

      double len = (p0->X-p1->X)*(p0->X-p1->X) + (p0->Y-p1->Y)*(p0->Y-p1->Y) +(p0->Z-p1->Z)*(p0->Z-p1->Z);
      if(len > lenla) {
        lenla = len;
        pedla = ped;
      }

    }               
    //printf("lenla %e %p : %d %d\n",lenla,pedla,EN_id((pEntity) E_vertex(pedla,0)),EN_id((pEntity) E_vertex(pedla,1)));
    EIter_delete(eit);  
    eit = M_edgeIter(mesh);  
    while ((ped = EIter_next(eit))) { 
      if(ped!=pedla) continue;
      /*les deux sommets de l'arete ped*/
      pVertex p0 = E_vertex(ped,0);
      pVertex p1 = E_vertex(ped,1);
      //pVertex pFix0,pFix1;
      std::vector<int> vectnod;
      int find = 0; 	
      /*sommets periodics ?*/
      void *temp_ptr0; 
      int isP0 = EN_getDataPtr((pEntity) p0 , tagPeriodic, &temp_ptr0);
      if(!isP0) continue;
      void *temp_ptr1; 
      int isP1 = EN_getDataPtr((pEntity) p1 , tagPeriodic, &temp_ptr1);
      if(!isP1) continue; 
      /*sommets deja tagges ?*/
      int move0;
      int isM0 = EN_getDataInt((pEntity) p0 ,tagMove, &move0);
      if(!(!isM0 || (move0 == -2))) continue;	  
      int move1;
      int isM1 = EN_getDataInt((pEntity) p1 ,tagMove, &move1);
      if(!(!isM1 || (move1 == -2))) continue;	  

      //printf("treat edge %d %d\n",EN_id((pEntity) p0),EN_id((pEntity) p1));    
      std::vector<std::pair<std::vector<int> , pVertex> > *recup0 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr0;
      std::vector<std::pair<std::vector<int> , pVertex> > *recup1 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;

      unsigned int j=0;
      for(j=0 ; j<(*recup0).size() ; j++) {
        std::vector<int> transfo0 = (*recup0)[j].first;
        pVertex pImg0 = (*recup0)[j].second;
        if(pImg0==p1) continue; 
        unsigned int j1;
        for(j1=0 ; j1<(*recup1).size() ; j1++) {
          std::vector<int> transfo1 = (*recup1)[j1].first;
          pVertex pImg1 = (*recup1)[j1].second;
          if(pImg1==p0) continue;  

          /*transfo0 == transfo1 ?*/
          assert(transfo1.size()==transfo0.size());
          unsigned int kt = 0;
          for(kt = 0 ; kt < transfo0.size(); kt++) {
            if(transfo0[kt] != transfo1[kt]) break;
          }
          if(kt!=transfo0.size()) continue; 
          /* pImg1-pImg0 existe ?*/
          if ( !E_exist(pImg0,pImg1) && !E_exist(pImg1,pImg0) ) continue;
          /* a-t-on deja trouve une arete img ?*/
          //printf("img edge %d %d\n",EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));    
          if(!find) {
            /*si le point img bouge, on bouge pas!*/
            int moveI0;
            int isMI0 = EN_getDataInt((pEntity) pImg0 ,tagMove, &moveI0);
            int moveI1;
            int isMI1 = EN_getDataInt((pEntity) pImg1 ,tagMove, &moveI1);
            if (isMI0 && (moveI0==1) ) continue;   	  
            if (isMI1 && (moveI1==1) ) continue;   	  
            find = 1;
            int move = 1;  
            EN_attachDataInt((pEntity) p0 ,tagMove, move);
            EN_attachDataInt((pEntity) p1 ,tagMove, move);
            move = -1;  
            EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
            EN_attachDataInt((pEntity) pImg1 ,tagMove, move);

            /*on garde les deux points img*/
            //pFix0 = pImg0;
            //pFix1 = pImg1;

            /*rajout de la transfo avec laquelle bouger les points*/
            for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
              vectnod.push_back(transfo0[kt]);
            }
            EN_attachDataPtr((pEntity) p0 , tagTransfo,new std::vector<int>(vectnod)); 
            EN_attachDataPtr((pEntity) p1 , tagTransfo,new std::vector<int>(vectnod)); 		      
            //printf("mark vertex : %d %d -- %d %d\n",EN_id((pEntity) p0),EN_id((pEntity) p1),EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));
            //printf("transfo : %d %d %d\n",transfo0[0],transfo0[1],transfo0[2]);
            break;
          } else {
            /*si le point bouge deja, on veut pas chger sa destination!!!*/
            int moveI0;
            int isMI0 = EN_getDataInt((pEntity) pImg0 ,tagMove, &moveI0);
            int moveI1;
            int isMI1 = EN_getDataInt((pEntity) pImg1 ,tagMove, &moveI1);
            if (isMI0 && (moveI0==1) ) continue;   	  
            if (isMI1 && (moveI1==1) ) continue;   	  

            int move = 1; 
            EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
            EN_attachDataInt((pEntity) pImg1 ,tagMove, move); 

            /*on sait que pImg0pImg1 doit aller sur pFix0pFix1 il faut trouver la transfo*/ 
            /*transfo = -transfo0 + vecnod  */
            std::vector<int> vecttrans;
            for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
              vecttrans.push_back(-transfo0[kt]);
            }
            for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
              vecttrans[kt] += vectnod[kt];
            }
            //printf("on applique la transfo %d %d %d\n",vecttrans[0],vecttrans[1],vecttrans[2]);
            //printf("mark vertex : %d %d \n",EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));
            EN_attachDataPtr((pEntity) pImg0 , tagTransfo,new std::vector<int>(vecttrans)); 
            EN_attachDataPtr((pEntity) pImg1 , tagTransfo,new std::vector<int>(vecttrans)); 		      
            break;			
          }	    
        } /*end j1*/
        if(j1==(*recup1).size()){
          //printf("arete not found\n");
          /* int move = -2;  
             EN_attachDataInt((pEntity) pImg0 ,tagMove, move);  */  	  
        }
      } /*end j*/	  
	//if(!find) {printf("on n'a pas trouve d'arete img donc on tag pas\n");       exit(0);}
    }  
    EIter_delete(eit);  
    eit = M_edgeIter(mesh);
    while ((ped = EIter_next(eit))) { 
      /*les deux sommets de l'arete ped*/
      pVertex p0 = E_vertex(ped,0);
      pVertex p1 = E_vertex(ped,1);
      //pVertex pFix0,pFix1;
      std::vector<int> vectnod;
      int find = 0;

      /*sommets periodics ?*/
      void *temp_ptr0; 
      int isP0 = EN_getDataPtr((pEntity) p0 , tagPeriodic, &temp_ptr0);
      if(!isP0) continue;
      void *temp_ptr1; 
      int isP1 = EN_getDataPtr((pEntity) p1 , tagPeriodic, &temp_ptr1);
      if(!isP1) continue; 
      /*sommets deja tagges ?*/
      int move0;
      int isM0 = EN_getDataInt((pEntity) p0 ,tagMove, &move0);
      if(!(!isM0 || (move0 == -2))) continue;	  
      int move1;
      int isM1 = EN_getDataInt((pEntity) p1 ,tagMove, &move1);
      if(!(!isM1 || (move1 == -2))) continue;	  

      //printf("treat edge %d %d\n",EN_id((pEntity) p0),EN_id((pEntity) p1));    
      std::vector<std::pair<std::vector<int> , pVertex> > *recup0 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr0;
      std::vector<std::pair<std::vector<int> , pVertex> > *recup1 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;
      unsigned int j=0;
      for(j=0 ; j<(*recup0).size() ; j++) {
        std::vector<int> transfo0 = (*recup0)[j].first;
        pVertex pImg0 = (*recup0)[j].second;
        if(pImg0==p1) continue; 
        unsigned int j1;
        for(j1=0 ; j1<(*recup1).size() ; j1++) {
          std::vector<int> transfo1 = (*recup1)[j1].first;
          pVertex pImg1 = (*recup1)[j1].second;
          if(pImg1==p0) continue;  

          /*transfo0 == transfo1 ?*/
          assert(transfo1.size()==transfo0.size());
          unsigned int kt = 0;
          for(kt = 0 ; kt < transfo0.size(); kt++) {
            if(transfo0[kt] != transfo1[kt]) break;
          }
          if(kt!=transfo0.size()) continue; 
          /* pImg1-pImg0 existe ?*/
          if ( !E_exist(pImg0,pImg1) && !E_exist(pImg1,pImg0) ) continue;
          /* a-t-on deja trouve une arete img ?*/
          //printf("img edge %d %d\n",EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));    
          if(!find) {
            /*si le point img bouge, on bouge pas!*/
            int moveI0;
            int isMI0 = EN_getDataInt((pEntity) pImg0 ,tagMove, &moveI0);
            int moveI1;
            int isMI1 = EN_getDataInt((pEntity) pImg1 ,tagMove, &moveI1);
            if (isMI0 && (moveI0==1) ) continue;   	  
            if (isMI1 && (moveI1==1) ) continue;   	  
            find = 1;
            int move = 1;  
            EN_attachDataInt((pEntity) p0 ,tagMove, move);
            EN_attachDataInt((pEntity) p1 ,tagMove, move);
            move = -1;  
            EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
            EN_attachDataInt((pEntity) pImg1 ,tagMove, move);

            /*on garde les deux points img*/
            //pFix0 = pImg0;
            //pFix1 = pImg1;

            /*rajout de la transfo avec laquelle bouger les points*/
            for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
              vectnod.push_back(transfo0[kt]);
            }
            EN_attachDataPtr((pEntity) p0 , tagTransfo,new std::vector<int>(vectnod)); 
            EN_attachDataPtr((pEntity) p1 , tagTransfo,new std::vector<int>(vectnod)); 		      
            //printf("mark vertex : %d %d -- %d %d\n",EN_id((pEntity) p0),EN_id((pEntity) p1),EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));
            break;
          } else {
            /*si le point bouge deja, on veut pas chger sa destination!!!*/
            int moveI0;
            int isMI0 = EN_getDataInt((pEntity) pImg0 ,tagMove, &moveI0);
            int moveI1;
            int isMI1 = EN_getDataInt((pEntity) pImg1 ,tagMove, &moveI1);
            if (isMI0 && (moveI0==1) ) continue;   	  
            if (isMI1 && (moveI1==1) ) continue;   	  
            int move = 1; 
            EN_attachDataInt((pEntity) pImg0 ,tagMove, move);
            EN_attachDataInt((pEntity) pImg1 ,tagMove, move); 

            /*on sait que pImg0pImg1 doit aller sur pFix0pFix1 il faut trouver la transfo*/ 
            /*transfo = -transfo0 + vecnod  */
            std::vector<int> vecttrans;
            for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
              vecttrans.push_back(-transfo0[kt]);
            }
            for(unsigned int kt = 0 ; kt < transfo0.size(); kt++) {
              vecttrans[kt] += vectnod[kt];
            }
            //printf("on applique la transfo %d %d %d\n",vecttrans[0],vecttrans[1],vecttrans[2]);
            //printf("mark vertex : %d %d \n",EN_id((pEntity) pImg0),EN_id((pEntity) pImg1));
            EN_attachDataPtr((pEntity) pImg0 , tagTransfo,new std::vector<int>(vecttrans)); 
            EN_attachDataPtr((pEntity) pImg1 , tagTransfo,new std::vector<int>(vecttrans)); 		      
            break;			
          }	    
        } /*end j1*/
        if(j1==(*recup1).size()){
          //printf("arete not found\n");
          /* int move = -2;  
             EN_attachDataInt((pEntity) pImg0 ,tagMove, move);  */  	  
        }
      } /*end j*/	  
	//if(!find) {printf("on n'a pas trouve d'arete img donc on tag pas\n");       exit(0);}
    }  
    EIter_delete(eit);
    return;
  }

  // -------------------------------------------------------------------
  void VertexTagMove_old(pMesh mesh,pMeshDataId tagMove,int imove) {
    pMeshDataId tagPeriodic     = MD_lookupMeshDataId("PeriodicPoint");

    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      void *temp_ptr; 
      int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
      if(isPeriodic) {
        /*****************************************************************************/
        /********nouvelle facon de bouger (sans pbs pour les tri-periodicite)*********/
        /*****************************************************************************/
        int move;
        int isMove = EN_getDataInt((pEntity) pv ,tagMove, &move);
        if((isMove && move==-1) || (isMove && move==1)) continue;	  
        /*****************************************************************************/
        /*****************************************************************************/
        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
        move = -1;
        unsigned int j=0;
        for(j=0 ; j<(*recup).size() ; j++) {
          std::vector<int> transfo = (*recup)[j].first;
          /*****************************************************************************/
          /********************ancienne facon*******************************************/
          /*****************************************************************************/		
          /*unsigned int k=0;
            for(k = 0 ; k<transfo.size() ; k++) {
            if(transfo[k]< 0) break;	
            }                    
            if(k==transfo.size()) move=1;*/  
          /*****************************************************************************/
		    
          /*on bouge que par rapport à x*/
	  /*  int somme = 0;
              for(int kk = 0 ; kk<transfo.size() ; kk++) {  
              //if(kk==imove) continue;
              somme +=  abs(transfo[kk] );
              }
              int inv = (imove>=0) ? 1 : -1;
              if(inv < 0) imove = abs(imove)-1; 
              for(k = imove ; k<imove+1 ; k++) {
              if((transfo[k] == 0) || (inv * transfo[k] < 0)) break;	
              }                 
              if(k==imove+1 && somme < 2) move=1;   */  
          /*****************************************************************************/
          /********nouvelle facon de bouger (sans pbs pour les tri-periodicite)*********/
          /*****************************************************************************/
          pVertex pImg = (*recup)[j].second;
          int movimg;
          int isM = EN_getDataInt((pEntity) pImg ,tagMove, &movimg);
          if(isM && movimg==1) continue;
    
          /*unsigned int k=0;
	    for(k = 0 ; k<transfo.size() ; k++) {
            if(transfo[k]< 0) break;	
            }                    
            if(k==transfo.size()) {
            move=1;
            break;
            } */
          move=1;
          break;
          /*****************************************************************************/
          /*****************************************************************************/
        }
        /*****************************************************************************/
        /********nouvelle facon de bouger (sans pbs pour les tri-periodicite)*********/
        /*****************************************************************************/
        if(move==1) {    
          //printf("on traite le point %d -> %d: j = %d / %d\n",EN_id((pEntity) pv),EN_id((pEntity) (*recup)[j].second),j,(*recup).size());
          unsigned int j2=0;
          for(j2=0 ; j2<(*recup).size() ; j2++) {
            std::vector<int> transfo = (*recup)[j2].first;
            pVertex pImg = (*recup)[j2].second;
            if(j2==j) {  
              int movimg;
              int isM = EN_getDataInt((pEntity) pImg ,tagMove, &movimg);    
              assert(!isM || (isM && movimg==-1) || (isM && movimg==-2));   
              movimg = -1;
              if(!isM)  EN_attachDataInt((pEntity) pImg ,tagMove, movimg);    
              continue;
            }
            /*le point pImg doit bouger sur le meme point que pv*/
            int movimg;
            int isM = EN_getDataInt((pEntity) pImg ,tagMove, &movimg); 
            if(!(!isM || (isM && movimg==-2) || (isM && movimg==1))) printf("point %d : %d %d\n",EN_id((pEntity) pImg),isM,movimg);   
            assert(!isM || (isM && movimg==-2) || (isM && movimg==1)); 
            movimg = 1;
            if(!isM)  EN_attachDataInt((pEntity) pImg ,tagMove, movimg);    
            // GCRMK: remove this part
            /****verifier que (*recup)[j].second existe dans les img de pImg*/
            pVertex pcheck = (*recup)[j].second;
            void *temp_ptrimg; 
            int isP = EN_getDataPtr((pEntity) pImg , tagPeriodic, &temp_ptrimg);
            assert(isP);
            std::vector<std::pair<std::vector<int> , pVertex> > *recupimg = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptrimg;
            unsigned int j3=0;
            for(j3=0 ; j3<(*recupimg).size() ; j3++) {
              pVertex pI2 = (*recup)[j3].second;
              if(pI2==pcheck) break;
            }
            assert(j3!=(*recupimg).size()); 
            /****end check*/
          }
        } else {
          move=-2;  
        }
        /*****************************************************************************/
        /*****************************************************************************/
        //printf("on traite le point %d : %d\n",EN_id((pEntity) pv),move);
        EN_attachDataInt((pEntity) pv ,tagMove, move);
      }
    }  
    VIter_delete(vit);  

    return;
  }

  // -------------------------------------------------------------------
  void VertexTagMoveInverse(pMesh mesh,pMeshDataId tagMove) {
    pMeshDataId tagPeriodic     = MD_lookupMeshDataId("PeriodicPoint");

    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      void *temp_ptr; 
      int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
      if(isPeriodic) {
        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
        int move = -1;
        for(unsigned int j=0 ; j<(*recup).size() ; j++) {
          std::vector<int> transfo = (*recup)[j].first;
          unsigned int k=0;
          for(k = 0 ; k<transfo.size() ; k++) {
            if(pv->X > 1) printf("%e %e %d\n",pv->X,pv->Y,transfo[k]);
            if(transfo[k] >= -1) break;
            //if((transfo[k] >= 0) || (transfo[k] >= -1 && transfo[k] <= 1)) break;	
            printf("tr %d\n",transfo[k]);
          }
          if(k==transfo.size()) move=1;
        }
        EN_attachDataInt((pEntity) pv ,tagMove, move);
      }
    }  
    VIter_delete(vit);  

    return;
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void MarkPeriodicTets(pMesh mesh,std::vector<std::vector<int> >& transfo,
                        pMeshDataId tagElt,pMeshDataId tagMove,pMeshDataId tagTransfo) {
  
    /* 1) define if periodic vertex can move : 1 : yes / -1 : no*/
    VertexTagMoveList(mesh,transfo,tagMove,tagTransfo);   
    /*ancienne version*/
    //VertexTagMove_old(mesh,tagMove,1);
  
    /* 2) marked tetras and vertex*/
    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      pVertex nod[4];
      pr->getNodes(nod);  
      int dest = 0;
      for(int i=0 ; i<4 ; i++) {
        int move; 
        int isPeriodic = EN_getDataInt((pEntity) nod[i] , tagMove, &move);  
        //GCRMK: remove this check
        if(isPeriodic && move==-2) printf("point %d argggg \n",EN_id((pEntity) nod[i])) ;
        if(isPeriodic) assert(move!=-2);
        if(isPeriodic && (move==1)) {
          dest = 1;	
        }
      }
      if(dest) EN_attachDataInt((pEntity) pr ,tagElt, dest);
    }  
    RIter_delete(rit);  

  }


  // -------------------------------------------------------------------
  int MarkGroupPeriodicTets(pMesh mesh,pMeshDataId tagElt,pMeshDataId tagMove) {
  
    int mark = 0;
  
    /* 1) define if periodic vertex can move : 1 : yes / -1 : no*/
    VertexTagMoveInverse(mesh,tagMove);
  
    /* 2) marked tetras and vertex*/
    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      pVertex nod[4];
      pr->getNodes(nod);  
      int dest = 0;
      for(int i=0 ; i<4 ; i++) {
        int move; 
        int isPeriodic = EN_getDataInt((pEntity) nod[i] , tagMove, &move);
        if(isPeriodic && (move==1)) {
 	  dest = 1;	
        }
      }
      if(dest) {
        EN_attachDataInt((pEntity) pr ,tagElt, dest);
        mark++;
      }
    }  
    RIter_delete(rit); 
    return mark; 

  }

  // -------------------------------------------------------------------
  void MarkPeriodicTriangles(pMesh mesh,std::vector<std::vector<int> >& transfo,pMeshDataId tagElt,pMeshDataId tagMove,pMeshDataId tagTransfo) {


    /* 1) define if periodic vertex can move : 1 : yes / -1 : no*/
    //VertexTagMove(mesh,tagMove,tagTransfo);
    VertexTagMoveList(mesh,transfo,tagMove,tagTransfo);
  
    /* 2) marked faces and vertex*/
    FIter fit = M_faceIter(mesh);
    pFace pface;  
    while ((pface = FIter_next(fit))) {
      pVertex nod[3];
      pface->getNodes(nod);  
      int dest = 0;
      for(int i=0 ; i<3 ; i++) {
        int move; 
        int isPeriodic = EN_getDataInt((pEntity) nod[i] , tagMove, &move);
        if(isPeriodic && (move==1)) {
 	  dest = 1;	
        }
      }
      if(dest) EN_attachDataInt((pEntity) pface ,tagElt, dest);
    }  
    FIter_delete(fit);
  }

  // -------------------------------------------------------------------

}

