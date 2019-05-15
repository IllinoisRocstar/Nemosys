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
#include "MeshDataBaseCommPeriodic.h"
#include "MeshDataBaseComm.h"
#include "MeshDataBaseParallelInterface.h"
#include "Mark.h"
#include "assert.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>

int ddebug = 0; 

namespace MAd {

  void MarkPeriodicEltVertex(pMesh mesh,pMeshDataId tagElt,pMeshDataId tagMove) {

    //Propagation du tagMove ie tagMove ==2
    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      int move; 
      int isMove = EN_getDataInt((pEntity) pv , tagMove,&move);
    
      if(!isMove) continue;
      if(move != 1) continue;
    
      int numed = V_numEdges(pv);
      for(int i=0 ; i<numed ; i++) {
        pEdge ped = V_edge(pv,i);
        pVertex pother = E_otherVertex(ped,pv);
        int movother; 
        int is = EN_getDataInt((pEntity) pother , tagMove,&movother);
      
        if(!is || (movother!=1)) EN_attachDataInt((pEntity) pother , tagMove,2);      
      } 
    }
    VIter_delete(vit);

    //tagMove == 10 : the vertex must be deleted
    //tagMove == -2 : the vertex don't move 
    vit = M_vertexIter(mesh);
    while ((pv = VIter_next(vit))) {
      int move; 
      int isMove = EN_getDataInt((pEntity) pv , tagMove,&move);
    
      if(!isMove) {
        EN_attachDataInt((pEntity) pv , tagMove,-2);   
      } else {
        if(abs(move)==1) continue;

        int numed = V_numEdges(pv);
        bool deleted = true;
        for(int i=0 ; i<numed ; i++) {
          pEdge ped = V_edge(pv,i);
          pVertex pother = E_otherVertex(ped,pv);
          int movother; 
          int is = EN_getDataInt((pEntity) pother , tagMove,&movother);
	
          if((!is) || (movother < 0)) {
            deleted = false;
          }
        }
        //if(deleted) EN_attachDataInt((pEntity) pv , tagMove,10);         
      }    
    }
    VIter_delete(vit);
  
  }
  bool compare(pEdge p1, pEdge p2);
  /*renvoie la bonne transfo pour les cas 3-periodic*/
  int EltMoveNew(const int dim,const pRegion pr,const pVertex* nod,pMeshDataId tagMove,pMeshDataId tagTransfo,std::vector<int> *vecttransfo) {
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
    int transformation = 0;
    //if(EN_id(R_vertex(pr,0))==5 || EN_id(R_vertex(pr,1))== 5 || EN_id(R_vertex(pr,2)) ==5 ||
    //							EN_id(R_vertex(pr,3))==5 )
    //printf("------ tet %d : %d %d %d %d\n",EN_id((pEntity) pr),EN_id(R_vertex(pr,0)),EN_id(R_vertex(pr,1)),EN_id(R_vertex(pr,2))
    //							,EN_id(R_vertex(pr,3)));  
    //if(EN_id(R_vertex(pr,0))==28 && EN_id(R_vertex(pr,1))== 111 && EN_id(R_vertex(pr,2)) ==4 &&
    //							EN_id(R_vertex(pr,3))==8 )
    //printf("****** tet %d : %d %d %d %d\n",EN_id((pEntity) pr),EN_id(R_vertex(pr,0)),EN_id(R_vertex(pr,1)),EN_id(R_vertex(pr,2))
    //													,EN_id(R_vertex(pr,3)));  
     
    /*on traite les aretes dans l'ordre decroissant (au cas ou plusieurs veulent bouger vers des dest diff)*/
    std::list<pEdge> listedge;
    for(int j=0 ; j<6 ; j++) {
      pEdge ped = R_edge(pr,j);   
      listedge.push_back(ped);
    } 
    //warning: anisotropic case!!
    listedge.sort(compare);

    /*boucle sur les aretes de pr*/
    for(std::list<pEdge>::iterator it = listedge.begin() ; it!=listedge.end() ; it++) {  
      pEdge ped = (*it);
    
      /*les deux sommets de l'arete ped*/
      pVertex p0 = E_vertex(ped,0);
      pVertex p1 = E_vertex(ped,1);
      if(EN_id(R_vertex(pr,0))==28 && EN_id(R_vertex(pr,1))== 111 && EN_id(R_vertex(pr,2)) ==4 &&
         EN_id(R_vertex(pr,3))==8 )
        printf("treat edge %d %d \n",EN_id((pEntity) p0),EN_id((pEntity) p1));    
    
      /*sommets periodics ?*/
      void *temp_ptr0; 
      int isP0 = EN_getDataPtr((pEntity) p0 , tagPeriodic, &temp_ptr0);
      if(!isP0) continue;
      void *temp_ptr1; 
      int isP1 = EN_getDataPtr((pEntity) p1 , tagPeriodic, &temp_ptr1);
      if(!isP1) continue;  
	
      void *tr0,*tr1; 
      int isT0 = EN_getDataPtr((pEntity) p0 ,tagTransfo, &tr0); 
      if(!isT0) continue;
      std::vector<int> *t0 = (std::vector<int> *) tr0;
      int isT1 = EN_getDataPtr((pEntity) p1 ,tagTransfo, &tr1);
      std::vector<int> *t1 = (std::vector<int> *) tr1;
      if(!isT1) continue;    
      unsigned int jj;
      for(jj=0 ; jj<(*t0).size() ; jj++) {   
        if((*t0)[jj] != (*t1)[jj]) break;   
      }
      if(jj!=(*t0).size()) {  
        //printf("pbs les 2 points ne bougent pas pareil!!! %d %d\n"
        //			,EN_id((pEntity) p0),EN_id((pEntity) p1));
        continue; 
      }
      //if(EN_id(R_vertex(pr,0))==28 && EN_id(R_vertex(pr,1))== 111 && EN_id(R_vertex(pr,2)) ==4 &&
      //							EN_id(R_vertex(pr,3))==8 )      
      //							printf("on bouge : %d %d %d\n",(*t0)[0],(*t0)[1],(*t0)[2]);    
      transformation = 1;
      //for(int kv=0 ; kv<transfo.size() ; kv++)
      (*vecttransfo) = *t0;
      return(transformation);
    } /*end j*/
    //printf("problemmmm\n");
    return transformation;
  } 
  int EltMove2dNew(const int dim,const pVertex* nod,pMeshDataId tagMove,pMeshDataId tagTransfo,std::vector<int> *vecttransfo) {
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
    int transformation = 0;

    for(int i=0 ; i<(dim+1) ; i++) {
      void *temp_ptr; 
      /*sommets periodics */
      int isPeriodic = EN_getDataPtr((pEntity) nod[i] , tagPeriodic, &temp_ptr);
      if(!isPeriodic) continue;
      if(ddebug) printf("test nod %d) \n",EN_id((pEntity) nod[i]));

      /*transfo*/                                               
      void *tr0; 
      int isT0 = EN_getDataPtr((pEntity) nod[i] ,tagTransfo, &tr0); 
      if(!isT0) continue;
      std::vector<int> *t0 = (std::vector<int> *) tr0;
      transformation = 1;
      //for(int kv=0 ; kv<transfo.size() ; kv++)
      (*vecttransfo) = *t0;
      return(transformation);

    }
    return transformation;
  }
  int EltMove(const int dim,const pVertex* nod,pMeshDataId tagMove,std::vector<int> *vecttransfo) {
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
    int transformation = 0;

 
    if(ddebug) printf("tetra %d %d %d %d\n",nod[0]->iD,nod[1]->iD,nod[2]->iD,nod[3]->iD);
    for(int i=0 ; i<(dim+1) ; i++) {
      void *temp_ptr; 
      int isPeriodic = EN_getDataPtr((pEntity) nod[i] , tagPeriodic, &temp_ptr);
      if(!isPeriodic) continue;
      if(ddebug) printf("test nod %d) \n",EN_id((pEntity) nod[i]));
      std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
      unsigned int j=0;
      for(j=0 ; j<(*recup).size() ; j++) {
        std::vector<int> transfo = (*recup)[j].first;
        /*****************************************************************************/
        /********************ancienne facon*******************************************/
        /*****************************************************************************/		
        unsigned int kk=0;
        for(kk = 0 ; kk<transfo.size() ; kk++) {
          if( transfo[kk]< 0) break;     
        }
        if(kk!=transfo.size()) continue;  
   	/*****************************************************************************/
        /*on bouge que par rapport à x*/   
        /**   int somme = 0;
              for(int kk = 0 ; kk<transfo.size() ; kk++) {
              //if(kk==imove) continue;
              somme +=  abs(transfo[kk] );
              } 
              int inv = (imove>=0) ? 1 : -1;    
              if(inv < 0) imove = abs(imove)-1; 
              for(kk = imove ; kk<imove+1 ; kk++) {
              if((transfo[kk] == 0) || (inv * transfo[kk] <0)) break;	
              }        
              if(kk!=imove+1 || somme >=2) continue;    */ 
        /*****************************************************************************/
        /********nouvelle facon de bouger (sans pbs pour les tri-periodicite)*********/ 
        /********on prend la seule transformation qui repond : ***********************/ 
        /******** tous les autres points img doivent etre en mvt *********************/
        /*****************************************************************************/

        /*****************************************************************************/
        /*****************************************************************************/

        if(ddebug) printf("on peut eventuellement le bouge par %d %d %d (point %d)\n",transfo[0],transfo[1],transfo[2],EN_id((pEntity) (*recup)[j].second));
        pVertex pImg = (*recup)[j].second;
        int move;
        int isMove = EN_getDataInt((pEntity) pImg ,tagMove, &move);
        assert(isMove);
        if(ddebug) printf("point img tagMove %d\n",move);      
        if(move==10 || move==1) continue;
        if(ddebug) printf("on teste de bouger sur le point %d (%d)\n",EN_id((pEntity) (*recup)[j].second),move);  
        /*******************************************************************************************************/
        /***** il faut verifier que s'il existe une img aux autres points par cette transfo elle soit fixe ****/
        /******************************************************************************************************/
        int k=0;//i+1;
        for(k=i+1 ; k<(dim+1) ; k++) { 
          if(k==i) continue;
          void *temp_ptr2; 
          int isPeriodic2 = EN_getDataPtr((pEntity) nod[k] , tagPeriodic, &temp_ptr2);
          if(!isPeriodic2) continue;
          if(ddebug) printf("on teste le point %d\n",EN_id((pEntity) nod[k]));
          std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
          unsigned int j2=0;
          for(j2=0 ; j2<(*recup2).size() ; j2++) {
            std::vector<int> transfo2 = (*recup2)[j2].first;
            assert(transfo2.size()==transfo.size());
            unsigned int k2=0;
            for(k2 = 0 ; k2<transfo2.size() ; k2++) { 
              if( transfo2[k2] != transfo[k2] ) break;     
            }
            if(k2!=transfo2.size()) continue;  
            pVertex pImg2 = (*recup2)[j2].second;
            int move2;
            int isMove2 = EN_getDataInt((pEntity) pImg2 ,tagMove, &move2);
            assert(isMove2);
            if(ddebug) printf("img %d (%d)\n",EN_id((pEntity) pImg2),move2);	  
            if(move2==10 || move2==1) break;	  
          }
          if(j2!=(*recup2).size()) break;//si on est arrive au bout de la boucle : ce noeud est ok on passe au suivant
        }
        if(k==(dim+1)) { //on peut bouger
          transformation = 1;
          //for(int kv=0 ; kv<transfo.size() ; kv++)
          (*vecttransfo) = transfo;
          break;
        }
      }  
      if(j!=(*recup).size()) break;//on peut bouger         
    }
    return transformation;
  }

  pFace findFace(pMesh mesh,const pVertex p1,const pVertex p2,const  pVertex p3) {
    pFace pface;
  
    pface =  F_exist(p1,p2,p3,0);
    if(!pface) pface = F_exist(p1,p2,p3,0);
    if(!pface) pface = F_exist(p1,p3,p2,0);
    if(!pface) pface = F_exist(p2,p1,p3,0);
    if(!pface) pface = F_exist(p2,p3,p1,0);
    if(!pface) pface = F_exist(p3,p1,p2,0);
    if(!pface) pface = F_exist(p3,p2,p1,0);
  
    assert(pface);
    return pface;

  }
  void TetraMove(pMesh mesh,pRegion pr,const pVertex* nodnew) {
    pGEntity pg = EN_whatIn(pr);
    int tag = GEN_tag(pg);
    int dim = GEN_type(pg); 
    if(dim==3) {
      pg = (pGEntity) GM_regionByTag(mesh->model,tag);      
    } else {
      printf("----pbs faces**** %d\n",dim);
    }

    pFace pface[4];
    pface[0] =  findFace(mesh,nodnew[0],nodnew[1],nodnew[2]);
    assert(pface[0]);
    pface[1] =  findFace(mesh,nodnew[0],nodnew[1],nodnew[3]);
    assert(pface[1]);
    pface[2] =  findFace(mesh,nodnew[1],nodnew[2],nodnew[3]);
    assert(pface[2]);
    pface[3] =  findFace(mesh,nodnew[0],nodnew[2],nodnew[3]);
    assert(pface[3]);
    //if((EN_id(nodnew[0])==5 || EN_id(nodnew[1])==5 ||EN_id(nodnew[2])==5 || EN_id(nodnew[3])==5) && 
    //		(EN_id(nodnew[0])==3 || EN_id(nodnew[1])==3 || EN_id(nodnew[2])==3) || EN_id(nodnew[3])==3)  
    //printf("new tet %d %d %d %d \n",EN_id(nodnew[0]),EN_id(nodnew[1]),EN_id(nodnew[2]),EN_id(nodnew[3]));
   
    pRegion prnew = M_createR(mesh,4,pface,pg);
    //printf("tetnew face : %p %p %p %p\n",prnew->f1,prnew->f2,prnew->f3,prnew->f4);		
    assert(prnew);
  
  
    return;
  }

  void TriangleMove(pMesh mesh,pFace pface,const pVertex* nodnew) {
    pGEntity pg = EN_whatIn(pface);
    int tag = GEN_tag(pg);
    int dim = GEN_type(pg); 
    if(dim==2) {
      pg = (pGEntity) GM_faceByTag(mesh->model,tag);
    } else if(dim==3) {
      pg = (pGEntity) GM_regionByTag(mesh->model,tag);      
    } else {
      printf("----pbs faces**** %d\n",dim);
    }
    pEdge pe[3];
    pe[0] =  E_exist(nodnew[0],nodnew[1]);
    assert(pe[0]);
    pe[1] =  E_exist(nodnew[1],nodnew[2]);
    assert(pe[1]);
    pe[2] =  E_exist(nodnew[0],nodnew[2]);
    //if(!pe[2]) pe[2] =  E_exist(nodnew[2],nodnew[0]);
    assert(pe[2]);
   
    pFace pfnew =F_exist(pe[0],pe[1],pe[2],0);
    if(!pfnew) pfnew = M_createF(mesh,3,pe,pg);
    //printf("tetnew face : %p\n",pfnew); 
   	
    assert(pfnew);
  

  
    return;
  }

  void TriangleDelete(pMesh mesh,const pVertex* nod) {
    int i1,i2,i3;
    for(int i=0 ; i<4 ; i++) {
      switch(i) {
      case 0 : 
        i1 = 0;
        i2 = 1;
        i3 = 2;
        break;
      case 1 : 
        i1 = 0;
        i2 = 1;
        i3 = 3;
        break;
      case 2 : 
        i1 = 1;
        i2 = 2;
        i3 = 3;
        break;
      case 3 : 
        i1 = 0;
        i2 = 2;
        i3 = 3;
        break;
      }
      pFace pface = F_exist(nod[i1],nod[i2],nod[i3],0);
      assert(pface);
      int num = F_numRegions(pface);
      //printf(" del? %d face %p \n",!num,pface);
      if (!num) M_removeFace(mesh,pface);
    }
    return;
  }

  void numP(const int cas,int* i1,int* i2) {
    switch(cas) {
    case 0 :
      *i1 = 0;
      *i2 = 1;
      break;  
    case 1 :
      *i1 = 1;
      *i2 = 2;
      break;  
    case 2 :
      *i1 = 2;
      *i2 = 0;
      break;  
    case 3 :
      *i1 = 0;
      *i2 = 3;
      break;  
    case 4 :
      *i1 = 1;
      *i2 = 3;
      break;  
    case 5 :
      *i1 = 2;
      *i2 = 3;
      break;  
    default : 
      throw;
      break;  
    }
    return;
  }

  void EdgeDelete(const int nbe,pMesh mesh,const pVertex* nod) {
    int i1,i2;
    for(int i = 0 ; i<nbe ; i++) {
      numP(i,&i1,&i2);
      pEdge pedge = E_exist(nod[i1],nod[i2]);
      //if(!pedge) pedge = E_exist(nod[i1],nod[i2]);
      assert(pedge);
      int num = E_numFaces(pedge); 
      if (!num) M_removeEdge(mesh,pedge);
    }
    return;
  }

  void EdgeMove(const int nbe,pMesh mesh,const pVertex* nod,const pVertex* nodnew) {
    int i1,i2;
    for(int i = 0 ; i<nbe ; i++) {
      numP(i,&i1,&i2);
   
      pEdge pedge = E_exist(nod[i1],nod[i2]);
      //if(!pedge) pedge = E_exist(nod[i2],nod[i1]);
      assert(pedge);
      pGEntity  pg = EN_whatIn(pedge);
      int dim  = GEN_type(pg);
      int tag  = GEN_tag(pg);
      if(dim==1) {
        pg = (pGEntity) GM_edgeByTag(mesh->model,tag);	   
      } else if(dim==2) {
        pg = (pGEntity) GM_faceByTag(mesh->model,tag);
      } else if(dim==3) {
        pg = (pGEntity) GM_regionByTag(mesh->model,tag);      
      } else {
        printf("----pbs**** %d\n",dim);
      }
      pEdge     e = E_exist(nodnew[i1],nodnew[i2]);  
      //if(EN_id(nod[i1])==4694 || EN_id(nod[i2])==4694 || EN_id(nodnew[i1])==4694 || EN_id(nod[i2])==4694)
      //	printf("edge %d  %d ---> %d %d\n",EN_id(nod[i1]),EN_id(nod[i2]),EN_id(nodnew[i1]),EN_id(nodnew[i2])) ;
      //printf("on cree edge : %d %d\n",EN_id((pEntity) nodnew[i1]),EN_id((pEntity) nodnew[i2]));
      //if(nodnew[i1]->iD==144 || nodnew[i2]->iD==144) printf("-------------------- on cree edge : %d %d\n",EN_id((pEntity) nodnew[i1]),EN_id((pEntity) nodnew[i2]));
      if (!e)   e = M_createE(mesh,nodnew[i1],nodnew[i2],pg);
    }
    return;
  }

  void VertexDelete(pMesh mesh, pMeshDataId tagMove, pMeshDataId tagTransfo){
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");

    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      int move;
      int isMove = EN_getDataInt((pEntity) pv ,tagMove, &move);
      assert(isMove);
      EN_deleteData((pEntity)pv, tagMove);
      //  if(!(move==1 || move==10)) continue;    
    
      //  printf("on delete le point %d ?\n",EN_id((pEntity) pv));
      int num = V_numEdges(pv);
      if(num){
        // if(move!=10) printf("point not deleted %d (%d)\n",EN_id((pEntity) pv),move);
        //assert(move==10);
        continue;
      }
      void *ttr; 
      int isT = EN_getDataPtr((pEntity) pv , tagTransfo, &ttr);
      if(isT) {
        std::vector<int> *rcup = (std::vector<int> *) ttr;
        delete rcup;
        EN_deleteData((pEntity) pv , tagTransfo);
      }
 
      void *temp_ptr; 
      int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
      if(isPeriodic) {
        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
        for(unsigned int j=0 ; j<(*recup).size() ; j++) {
          pVertex ptemp = (*recup)[j].second;
          void *temp_ptr2; 
          int is = EN_getDataPtr((pEntity) ptemp , tagPeriodic, &temp_ptr2);
          assert(is);
          std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
          if((*recup2).size()==1) {
            EN_deleteData((pEntity) ptemp , tagPeriodic);
          } else {
            std::vector<std::pair<std::vector<int> , pVertex> > vectnod;
            for(unsigned int k=0 ; k<(*recup2).size() ; k++) {
              if((*recup2)[k].second != pv) vectnod.push_back(std::pair<std::vector<int> , pVertex>((*recup2)[k].first,(*recup2)[k].second));
            }	
            EN_attachDataPtr((pEntity) ptemp , tagPeriodic, 
                             new std::vector<std::pair<std::vector<int> , pVertex> >(vectnod)); 
          }
        }      
        delete recup;
        EN_deleteData((pEntity) pv , tagPeriodic);
      }
      
      M_removeVertex(mesh,pv);
    
    }
    VIter_delete(vit);

    return;
  }

  void VertexMove(const int nbv,pMesh mesh,const pVertex* nod,const std::vector<int>& vecttransfo,pMeshDataId tagMove,
                  pVertex* nodnew,MDB_DataExchangerPeriodic &deperiodic) {
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
      
    std::vector<int> invvecttransfo;  
    for(unsigned k=0 ; k<vecttransfo.size() ; k++) {
      invvecttransfo.push_back(-vecttransfo[k]);
    }
  
    for(int i=0 ; i<nbv ; i++) {
      void *temp_ptr; 
      int isPeriodic = EN_getDataPtr((pEntity) nod[i] , tagPeriodic, &temp_ptr);
      if(!isPeriodic) {
        //creation du point, nod[i] et nodnew[i] devient periodic
        pGEntity  pg  = EN_whatIn(nod[i]);
        pGEntity pent;
        int dim  = GEN_type(pg);
        int tag  = GEN_tag(pg);
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
        double X,Y,Z;
        X = nod[i]->X;
        Y = nod[i]->Y;
        Z = nod[i]->Z;
        for(int k = 0 ; k<deperiodic.nbRefPeriodic() ; k++) {
          if(vecttransfo[k]){
            int inv1 = (vecttransfo[k] < 0) ? 1 : 0;
            for(int nb = 0 ; nb<abs(vecttransfo[k]) ; nb++) {
              deperiodic.fperiodic(inv1,X,Y,Z,k+1,&X,&Y,&Z);
            }  
          }
        }
        nodnew[i] = M_createV(mesh,X,Y,Z,-1,pent);
        if(ddebug) {
          printf("nod vient de %d avec %d %d\n",EN_id((pEntity) nod[i]),vecttransfo[0],vecttransfo[1]);
          printf("coor nod %e %e\n",nod[i]->X,nod[i]->Y);
          printf("coor nodnew %e %e\n",nodnew[i]->X,nodnew[i]->Y);	
        }
        std::vector<std::pair<std::vector<int> , pVertex> > vectnodperiodic;
        vectnodperiodic.push_back(std::pair<std::vector<int> , pVertex>(vecttransfo,nodnew[i]));
        EN_attachDataPtr((pEntity) nod[i] , tagPeriodic, 
                         new std::vector<std::pair<std::vector<int> , pVertex> >(vectnodperiodic)); 	
        std::vector<std::pair<std::vector<int> , pVertex> > vectnodnewperiodic;

        vectnodnewperiodic.push_back(std::pair<std::vector<int> , pVertex>(invvecttransfo,nod[i]));
        EN_attachDataPtr((pEntity) nodnew[i] , tagPeriodic, 
                         new std::vector<std::pair<std::vector<int> , pVertex> >(vectnodnewperiodic)); 
        EN_attachDataInt((pEntity) nodnew[i] , tagMove,2); 
        if(ddebug) printf("creation non periodique %d -> %d (2)\n",EN_id((pEntity) nod[i]),EN_id((pEntity) nodnew[i]));					
 
      } else {
        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
        unsigned int j=0;
        for(j=0 ; j<(*recup).size() ; j++) {
          std::vector<int> transfo = (*recup)[j].first;
          assert(transfo.size()==vecttransfo.size());
          unsigned int kt = 0;
          for(kt = 0 ; kt < vecttransfo.size(); kt++) {
            if(vecttransfo[kt] != transfo[kt]) break;
          }
          if(kt!=vecttransfo.size()) continue;  

          if((EN_id(nod[0])==28 && EN_id(nod[1])==111 && EN_id(nod[2])==4 && EN_id(nod[3])==8))  
            printf("le correspondant de %d est %d par %d %d\n",EN_id((pEntity) nod[i]),EN_id((pEntity)(*recup)[j].second),
                   vecttransfo[0],vecttransfo[1]);
          if(ddebug) printf("le correspondant de %d est %d par %d %d\n",EN_id((pEntity) nod[i]),EN_id((pEntity)(*recup)[j].second),
                            vecttransfo[0],vecttransfo[1]);	
          //le point existe
          nodnew[i] = (*recup)[j].second;
          break;
        }
        if(j==(*recup).size()) {
          //printf("ARGGGGGGGGGGGGGG point %d (%e %e) par %d %d\n",EN_id((pEntity) nod[i]),nod[i]->X,nod[i]->Y,
          //					vecttransfo[0],vecttransfo[1]);    
          //le point n'existe pas mais nod est periodic
          pGEntity  pg  = EN_whatIn(nod[i]);
          pGEntity pent;
          int dim  = GEN_type(pg);
          int tag  = GEN_tag(pg);
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
          double X,Y,Z;
          X = nod[i]->X;
          Y = nod[i]->Y;
          Z = nod[i]->Z;
          for(int k = 0 ; k<deperiodic.nbRefPeriodic() ; k++) {
            if(vecttransfo[k]){
              int inv1 = (vecttransfo[k] < 0) ? 1 : 0;
              for(int nb = 0 ; nb<abs(vecttransfo[k]) ; nb++) {  
                deperiodic.fperiodic(inv1,X,Y,Z,k+1,&X,&Y,&Z);
              }	  
            }
          }
          nodnew[i] = M_createV(mesh,X,Y,Z,-1,pent);
          //printf("on cree le point %p : %e %e \n",nodnew[i],X,Y);
          std::vector<std::pair<std::vector<int> , pVertex> > vectnodnewperiodic;
          //printf("1) on ajoute le point %d (%e %e) a la liste de %d (%e %e) par %d %d\n",EN_id((pEntity) nod[i])
          //			,nod[i]->X,nod[i]->Y
          //			,EN_id((pEntity) nodnew[i]),nodnew[i]->X,nodnew[i]->Y
          //			,invvecttransfo[0],invvecttransfo[1]);
          vectnodnewperiodic.push_back(std::pair<std::vector<int> , pVertex>(invvecttransfo,nod[i]));
          for(j=0 ; j<(*recup).size() ; j++) {
            std::vector<int> transfo = (*recup)[j].first;
            pVertex ptemp = (*recup)[j].second;
            void *temp_ptr2; 
            int is = EN_getDataPtr((pEntity) ptemp , tagPeriodic, &temp_ptr2); 
            assert(is);
            //#warning: calcule la bonne transformation	  
            std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
#ifdef DEBUG
            unsigned int k=0;
            for(k=0 ; k<(*recup2).size() ; k++) {
              std::vector<int> transfo2 = (*recup2)[k].first; 
              int kk=0;                                    
              for(kk = 0 ; kk<deperiodic.nbRefPeriodic() ; kk++) {       
                if(transfo[kk]!=(-1) * transfo2[kk]) break;
              }
              if(kk<deperiodic.nbRefPeriodic()) continue;  
              //printf(" %p (%e %e %e) -- %p (%e %e %e)\n",(*recup2)[k].second,(*recup2)[k].second->X,(*recup2)[k].second->Y,(*recup2)[k].second->Z,nod[i],nod[i]->X,nod[i]->Y,nod[i]->Z);
              assert((*recup2)[k].second == nod[i]); 
              break;
            }        
            assert(k!=(*recup2).size());
#endif
            std::vector<int> newtransfo;
            std::vector<int> invnewtransfo;
            for(unsigned int kn = 0 ; kn < transfo.size() ; kn++){
              invnewtransfo.push_back(vecttransfo[kn] - transfo[kn]);
              newtransfo.push_back(-(vecttransfo[kn] - transfo[kn]));
            }
            //printf("2) on ajoute le point %d (%e %e) a la liste de %d (%e %e) par %d %d\n",EN_id((pEntity) ptemp)
            //   			,ptemp->X,ptemp->Y
            //   			,EN_id((pEntity) nodnew[i]),nodnew[i]->X,nodnew[i]->Y
            //   			,newtransfo[0],newtransfo[1]);
            vectnodnewperiodic.push_back(std::pair<std::vector<int> , pVertex>(newtransfo,ptemp));
            //printf("3) on ajoute le point %d (%e %e) a la liste de %d (%e %e) par %d %d\n",EN_id((pEntity) nodnew[i])
            //  			,nodnew[i]->X,nodnew[i]->Y
            //  			,EN_id((pEntity) ptemp),ptemp->X,ptemp->Y
            //  			,invnewtransfo[0],invnewtransfo[1]);
            (*recup2).push_back(std::pair<std::vector<int> , pVertex>(invnewtransfo,nodnew[i]));
          }  
          EN_attachDataPtr((pEntity) nodnew[i] , tagPeriodic, 
                           new std::vector<std::pair<std::vector<int> , pVertex> >(vectnodnewperiodic)); 
          EN_attachDataInt((pEntity) nodnew[i] , tagMove,2); 
          (*recup).push_back(std::pair<std::vector<int> , pVertex>(vecttransfo,nodnew[i]));
          //printf("4) on ajoute le point %d (%e %e) a la liste de %d (%e %e) par %d %d\n",EN_id((pEntity) nodnew[i])
          //			,nodnew[i]->X,nodnew[i]->Y
          //			,EN_id((pEntity) nod[i]),nod[i]->X,nod[i]->Y
          //			,vecttransfo[0],vecttransfo[1]);
							     
        } else {
          //le point existe rien a faire pour l'instant
        }
      }//end if periodic
    }// end for i
    return;
  }

  void MovePeriodicTriangles(pMesh mesh,pMeshDataId tagElt,pMeshDataId tagMove,pMeshDataId tagTransfo, MDB_DataExchanger &de,
                             MDB_DataExchangerPeriodic &deperiodic) {

    int nmodif = 0;
    EIter eit;
    pEdge pedge;  
  
    do {
      printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ while : nb modif = %d\n",nmodif);
      //if(nmodif) break;
      nmodif = 0;
      eit = M_edgeIter(mesh);
      while ((pedge = EIter_next(eit))) {
        pVertex p[2];
        p[0] = pedge->p1;
        p[1] = pedge->p2;
        int move1,move2;
        int isMove1 = EN_getDataInt((pEntity) p[0] ,tagMove, &move1);
        int isMove2 = EN_getDataInt((pEntity) p[1] ,tagMove, &move2);
        assert(isMove1 && isMove2);
        if(!(move1 == 1 || move1 == 10 || move2==1 || move2 == 10)) continue;
        //l'arete doit bouger
        if(E_numFaces(pedge)!=1) continue;
        pFace pface = E_face(pedge,0);
        assert(pface);
        pVertex nod[3];
        pface->getNodes(nod);

        int tmp; 
        int isChange = EN_getDataInt((pEntity) pface , tagElt, &tmp);
        if(!isChange) continue;

        //       printf("test edge ------------- : %d (%d) -- %d (%d)\n",EN_id((pEntity)p[0]),move1,EN_id((pEntity)p[1]),move2);
        //       printf("tr(%d %d %d) \n\n",EN_id((pEntity)nod[0]),EN_id((pEntity)nod[1]),EN_id((pEntity)nod[2]));    
     
        //test si on peut bouger le triangle issu de cette arete
        std::vector<int> vecttransfo;
        int transformation = EltMove2dNew(2,nod,tagMove,tagTransfo,&vecttransfo);
        if(!transformation) continue; 
        //printf("------------- edge : %d (%d) -- %d (%d)\n",EN_id((pEntity)p[0]),move1,EN_id((pEntity)p[1]),move2);
        //printf("on a le droit de bouger (%d %d %d) selon %d %d\n\n",EN_id((pEntity)nod[0]),EN_id((pEntity)nod[1]),EN_id((pEntity)nod[2]),
        //					vecttransfo[0],vecttransfo[1]);    
      
       
        //si oui on le bouge => creation des points, edges, face + delete face, edges, points
        // + mise a jour des tagPeriodic et des tagMove (ie points deviennent fixes : -2)
     
        //bouge des points : s'il n'existe pas on le cree, sinon on donne son pointeur
        pVertex nodnew[3];
        VertexMove(3,mesh,nod,vecttransfo,tagMove,nodnew,deperiodic);
     
        //creation des edges
        EdgeMove(3,mesh,nod,nodnew);
        //printf("on bouge les points (%f %f) (%f %f) (%f %f) \n",nod[0]->X,nod[0]->Y,
        //   							nod[1]->X,nod[1]->Y,   
        //   							nod[2]->X,nod[2]->Y);
        //printf("vers     les points (%f %f) (%f %f) (%f %f) \n",nodnew[0]->X,nodnew[0]->Y,
        //   							nodnew[1]->X,nodnew[1]->Y,   
        //   							nodnew[2]->X,nodnew[2]->Y);
        //creation/delete des triangles
        TriangleMove(mesh,pface,nodnew);
        int is = EN_getDataInt((pEntity) pface , tagElt, &tmp);
        assert(is);	 
        EN_deleteData((pEntity) pface , tagElt);
        M_removeFace(mesh,pface);  
   
        EdgeDelete(3,mesh,nod); 
      
        nmodif++;
      }
      EIter_delete(eit);
    } while (nmodif);
  
    //check si plus de tr a bouger

  }

  void MovePeriodicTetras(pMesh mesh,pMeshDataId tagElt,pMeshDataId tagMove ,pMeshDataId tagTransfo,
                          MDB_DataExchanger &de,
                          MDB_DataExchangerPeriodic &deperiodic) {

    int nmodif = 0;
    FIter fit;
    pFace pface;  
  
    do {
      printf(" Tets $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ while : nb modif = %d\n",nmodif);
      //M_writeSMS(mesh,"chk.msh",10);
      //CheckMshPeriodic(mesh);
    
      //     //if(nmodif) break;
      nmodif = 0;
  
      //   //check adj
      //     fit = M_faceIter(mesh);
      //   while ((pface = FIter_next(fit))) {
      //     assert(!pface->deleted);
      //     int num = F_numRegions(pface);
      //     if(!num) {
      //       printf("face %p dim %d\n",pface,GEN_type(pface->g));
      //       exit(0);
      //       continue;
      //     }
      //   }
      //   FIter_delete(fit);
  
  
      fit = M_faceIter(mesh);
      while ((pface = FIter_next(fit))) {
       
        if(pface->deleted) continue;
        pVertex p[3];
        pface->getNodes(p);

        int move1,move2,move3;
        int isMove1 = EN_getDataInt((pEntity) p[0] ,tagMove, &move1);
        int isMove2 = EN_getDataInt((pEntity) p[1] ,tagMove, &move2);
        int isMove3 = EN_getDataInt((pEntity) p[2] ,tagMove, &move3);
        assert(isMove1 && isMove2 && isMove3);
        //	if((p[0]->iD==5 || p[1]->iD==5 || p[2]->iD==5) && (p[0]->iD==3 || p[2]->iD==3 || p[1]->iD==3))  
        //		printf("*** face %d %d %d : %d %d %d\n",p[0]->iD,p[1]->iD,p[2]->iD,move1,move2,move3) ;
	
        if(!(move1 == 1 || move1 == 10 || move2==1 || move2 == 10 || move3==1 || move3 == 10 )) continue;
        if(F_numRegions(pface)<1) printf("face %p pbs\n",(void*)pface);  
        assert(F_numRegions(pface)>=1);
        if(F_numRegions(pface)==2) continue;

        if((p[0]->iD==0 && p[1]->iD==0 && p[2]->iD==0) ||
           (p[0]->iD==0 && p[2]->iD==0 && p[1]->iD==0) ||
           (p[1]->iD==0 && p[2]->iD==0 && p[0]->iD==0) ||
           (p[1]->iD==0 && p[0]->iD==0 && p[2]->iD==0) ||
           (p[2]->iD==0 && p[1]->iD==0 && p[0]->iD==0) ||
           (p[2]->iD==0 && p[0]->iD==0 && p[1]->iD==0)) ddebug = 1;
        else if(p[0]->iD==0 || p[1]->iD==0 || p[2]->iD==0)ddebug = 1;
        else ddebug = 0;  
 
        pRegion pr = F_region(pface,0);
        assert(pr);
        pVertex nod[4];
        pr->getNodes(nod);

        int tmp; 
        int isChange = EN_getDataInt((pEntity) pr , tagElt, &tmp);
        if(!isChange) continue; 

        if(ddebug) printf("test face ------------- : %d (%d) -- %d (%d) -- %d(%d)\n",EN_id((pEntity)p[0]),move1,
                          EN_id((pEntity)p[1]),move2,
                          EN_id((pEntity)p[2]),move3);
     
        if(ddebug) printf("tet : %d %d %d %d\n",EN_id((pEntity)nod[0]),EN_id((pEntity)nod[1]),
                          EN_id((pEntity)nod[2]),EN_id((pEntity)nod[3]));
        if(ddebug) printf("tet face : %p %p %p %p\n",(void*)(pr->getFace(0)),(void*)(pr->getFace(1)),(void*)(pr->getFace(2)),(void*)(pr->getFace(3)));		
		
        //test si on peut bouger le tetra issu de cette face
        std::vector<int> vecttransfo;
        int transformation = EltMoveNew(3,pr,nod,tagMove,tagTransfo,&vecttransfo); 
        /*ancienne version*/ 
        //int transformation = EltMove(3,nod,tagMove,&vecttransfo);

        if(ddebug) printf("on bouge ? %d\n",transformation);
        if(!transformation) continue;       
       
      
        //Orientation!!!!!!
       
        pVertex nodnew[4];
        VertexMove(4,mesh,nod,vecttransfo,tagMove,nodnew,deperiodic);
	if((EN_id(nodnew[0])==23 && EN_id(nodnew[1])==148 && EN_id(nodnew[2])==3 && EN_id(nodnew[3])==5))  
          printf("tet %d %d %d %d: %d?\n",EN_id(nod[0]),EN_id(nod[1]),EN_id(nod[2]),EN_id(nod[3]),isChange);
        // printf("tetnew : %d %d %d %d\n",EN_id((pEntity)nodnew[0]),EN_id((pEntity)nodnew[1]),
        //       		EN_id((pEntity)nodnew[2]),EN_id((pEntity)nodnew[3]));
        //      
        //creation des edges
        EdgeMove(6,mesh,nod,nodnew);
     
        //creation/delete des triangles
        for(int i=0 ; i<4 ; i++) {
          pVertex nodtmp[3];
          switch(i) {
          case 0 : 
            nodtmp[0] = nodnew[0];
            nodtmp[1] = nodnew[1];
            nodtmp[2] = nodnew[2];
            break;
          case 1 : 
            nodtmp[0] = nodnew[0];
            nodtmp[1] = nodnew[1];
            nodtmp[2] = nodnew[3];
            break;
          case 2 : 
            nodtmp[0] = nodnew[1];
            nodtmp[1] = nodnew[2];
            nodtmp[2] = nodnew[3];
            break;
          case 3 : 
            nodtmp[0] = nodnew[0];
            nodtmp[1] = nodnew[2];
            nodtmp[2] = nodnew[3];
            break;
          }
          TriangleMove(mesh,pface,nodtmp);
        }

        TetraMove(mesh,pr,nodnew); 
      
        int is = EN_getDataInt((pEntity) pr , tagElt, &tmp);
        assert(is);	 
        EN_deleteData((pEntity) pr , tagElt);
        M_removeRegion(mesh,pr);  
     
        TriangleDelete(mesh,nod);
        EdgeDelete(6,mesh,nod); 
      
      
      
        nmodif++;
      }
      FIter_delete(fit);

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$MAJ des refs
      fit = M_faceIter(mesh);
      while ((pface = FIter_next(fit))) {
        assert(!pface->deleted);
        int num = F_numRegions(pface);
        if(!num) {
          printf("face %p dim %d\n",(void*)pface,GEN_type(pface->g));
          M_removeFace(mesh,pface); 
          continue;
        }
        int dim = GEN_type(pface->g); 
        if(dim == 2) {
          int num = F_numRegions(pface);
          if(num == 2) {
            //puts("enlever la ref");
            pRegion pr = F_region(pface,0);
            //	pGEntity pg  = EN_whatIn(pface);
            pGEntity pgr = EN_whatIn(pr);
            int tagr = GEN_tag(pgr);
            int dimr = GEN_type(pgr); 
            if(dimr==3) {
              pface->g = (pGEntity) GM_regionByTag(mesh->model,tagr);   
              assert(   GEN_type(pface->g)==3);
            } else {
              printf("----pbs faces**** %d\n",dim);
            }

          } else if(!num) {
            M_removeFace(mesh,pface); 
          }
          if(!(F_numRegions(pface)==1 || GEN_type(pface->g)!=2)) {
            printf("face %p %d %d %d\n",(void*)pface,F_numRegions(pface),GEN_type(pface->g),pface->deleted);
          }
          assert(F_numRegions(pface)==1 || GEN_type(pface->g)!=2 || pface->deleted);
        }
      }
      FIter_delete(fit);
      EIter eit = M_edgeIter(mesh);
      pEdge ped;
      while ((ped = EIter_next(eit))) {
        int num = E_numFaces(ped);
        if(!num) {
          printf("edge dim %d\n",GEN_type(ped->g));
          M_removeEdge(mesh,ped); 
          continue;
        }

      }
      EIter_delete(eit);
      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ End MAJ des refs


    } while (nmodif);
  
  }

  void PeriodicInterfaceMigration(pMesh mesh,pMeshDataId tagElt,pMeshDataId tagMove, 
                                  pMeshDataId tagTransfo,MDB_DataExchanger &de,
                                  MDB_DataExchangerPeriodic &deperiodic) {

    int dim = (mesh->tets.empty()) ? 2 : 3;

    /* 1) marked vertex*/
    MarkPeriodicEltVertex(mesh,tagElt,tagMove);

    /* 2) move elements*/
    if(dim==2) MovePeriodicTriangles(mesh,tagElt,tagMove,tagTransfo,de,deperiodic);
    else  MovePeriodicTetras(mesh,tagElt,tagMove,tagTransfo,de,deperiodic);
  
    /* 3) delete vertex*/
    VertexDelete(mesh,tagMove,tagTransfo);
  
  }
 
 
 
 
 
 
 
 
 
 
 
 
 
  int EltMoveInverse(const int dim,const pVertex* nod,pMeshDataId tagMove,std::vector<int> *vecttransfo) {
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
    int transformation = 0;
  
    for(int i=0 ; i<(dim+1) ; i++) {
      void *temp_ptr; 
      int isPeriodic = EN_getDataPtr((pEntity) nod[i] , tagPeriodic, &temp_ptr);
      if(!isPeriodic) continue;
      //printf("test nod %d)\n",EN_id((pEntity) nod[i]));
      std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
      unsigned int j=0;
      for(j=0 ; j<(*recup).size() ; j++) {
        std::vector<int> transfo = (*recup)[j].first;
        unsigned int kk=0;
        for(kk = 0 ; kk<transfo.size() ; kk++) {
          if( transfo[kk]> 0) break;     
        }
        if(kk!=transfo.size()) continue;
        //printf("on peut eventuellement le bouge par %d %d (point %d)\n",transfo[0],transfo[1],EN_id((pEntity) (*recup)[j].second));
        pVertex pImg = (*recup)[j].second;
        int move;
        int isMove = EN_getDataInt((pEntity) pImg ,tagMove, &move);
        assert(isMove);
        if(move==10 || move==1) continue;
        //printf("on teste de bouger sur le point %d (%d)\n",EN_id((pEntity) (*recup)[j].second),move);
        int k=i+1;
        for(k=i+1 ; k<(dim+1) ; k++) {
          void *temp_ptr2; 
          int isPeriodic2 = EN_getDataPtr((pEntity) nod[k] , tagPeriodic, &temp_ptr2);
          if(!isPeriodic2) continue;
          //printf("on teste le point %d\n",EN_id((pEntity) nod[k]));
          std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
          unsigned int j2=0;
          for(j2=0 ; j2<(*recup2).size() ; j2++) {
            std::vector<int> transfo2 = (*recup2)[j2].first;
            assert(transfo2.size()==transfo.size());
            unsigned int k2=0;
            for(k2 = 0 ; k2<transfo2.size() ; k2++) {
              if( transfo2[k2] != transfo[k2]) break;     
            }
            if(k2!=transfo2.size()) continue;
            pVertex pImg2 = (*recup2)[j2].second;
            int move2;
            int isMove2 = EN_getDataInt((pEntity) pImg2 ,tagMove, &move2);
            assert(isMove2);
            //printf("img %d (%d)\n",EN_id((pEntity) pImg2),move2);	  
            if(move2==10 || move2==1) break;	  
          }
          if(j2!=(*recup2).size()) break;//si on est arrive au bout de la boucle : ce noeud est ok on passe au suivant
        }
        if(k==(dim+1)) { //on peut bouger
          transformation = 1;
          //for(int kv=0 ; kv<transfo.size() ; kv++)
          (*vecttransfo) = transfo;
          break;
        }
      }  
      if(j!=(*recup).size()) break;//on peut bouger         
    }
    return transformation;
  } 
 
 
 
 
 
  void MovePeriodicTetrasInverse(pMesh mesh,pMeshDataId tagElt,pMeshDataId tagMove, MDB_DataExchanger &de,
                                 MDB_DataExchangerPeriodic &deperiodic) {

    int nmodif = 0;
    FIter fit;
    pFace pface;  
  
    do {
      printf("Inverse $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ while : nb modif = %d\n",nmodif);
      //     //if(nmodif) break;
      nmodif = 0;
  
      fit = M_faceIter(mesh);
      while ((pface = FIter_next(fit))) {
       
        if(pface->deleted) continue;
        pVertex p[3];
        pface->getNodes(p);
        int move1,move2,move3;
        int isMove1 = EN_getDataInt((pEntity) p[0] ,tagMove, &move1);
        int isMove2 = EN_getDataInt((pEntity) p[1] ,tagMove, &move2);
        int isMove3 = EN_getDataInt((pEntity) p[2] ,tagMove, &move3);
        assert(isMove1 && isMove2 && isMove3);
        if(!(move1 == 1 || move1 == 10 || move2==1 || move2 == 10 || move3==1 || move3 == 10 )) continue;
        if(F_numRegions(pface)==2) continue;
       
        pRegion pr = F_region(pface,0);
        assert(pr);
        pVertex nod[4];
        pr->getNodes(nod);

        int tmp; 
        int isChange = EN_getDataInt((pEntity) pr , tagElt, &tmp);
        if(!isChange) continue; 

        //printf("test face ------------- : %d (%d) -- %d (%d) -- %d(%d)\n",EN_id((pEntity)p[0]),move1,
        //				EN_id((pEntity)p[1]),move2,
        //				EN_id((pEntity)p[2]),move3);
     
        //printf("tet : %d %d %d %d\n",EN_id((pEntity)nod[0]),EN_id((pEntity)nod[1]),
        //      		EN_id((pEntity)nod[2]),EN_id((pEntity)nod[3]));
        //printf("tet face : %p %p %p %p\n",pr->f1,pr->f2,pr->f3,pr->f4);		
		
        //test si on peut bouger le tetra issu de cette face
        std::vector<int> vecttransfo;
        int transformation = EltMoveInverse(3,nod,tagMove,&vecttransfo);
        //printf("on bouge ? %d\n",transformation);
        if(!transformation) continue;       
       
       
        //Orientation!!!!!!
       
        pVertex nodnew[4];
        VertexMove(4,mesh,nod,vecttransfo,tagMove,nodnew,deperiodic);
        // printf("tetnew : %d %d %d %d\n",EN_id((pEntity)nodnew[0]),EN_id((pEntity)nodnew[1]),
        //       		EN_id((pEntity)nodnew[2]),EN_id((pEntity)nodnew[3]));
        //      
        //creation des edges
        EdgeMove(6,mesh,nod,nodnew);
     
        //creation/delete des triangles
        for(int i=0 ; i<4 ; i++) {
          pVertex nodtmp[3];
          switch(i) {
          case 0 : 
            nodtmp[0] = nodnew[0];
            nodtmp[1] = nodnew[1];
            nodtmp[2] = nodnew[2];
            break;
          case 1 : 
            nodtmp[0] = nodnew[0];
            nodtmp[1] = nodnew[1];
            nodtmp[2] = nodnew[3];
            break;
          case 2 : 
            nodtmp[0] = nodnew[1];
            nodtmp[1] = nodnew[2];
            nodtmp[2] = nodnew[3];
            break;
          case 3 : 
            nodtmp[0] = nodnew[0];
            nodtmp[1] = nodnew[2];
            nodtmp[2] = nodnew[3];
            break;
          }
          TriangleMove(mesh,pface,nodtmp);
        }

        TetraMove(mesh,pr,nodnew); 
      
        int is = EN_getDataInt((pEntity) pr , tagElt, &tmp);
        assert(is);	 
        EN_deleteData((pEntity) pr , tagElt);
        M_removeRegion(mesh,pr);  
     
        TriangleDelete(mesh,nod);
        EdgeDelete(6,mesh,nod); 
      
        nmodif++;
      }
      FIter_delete(fit);
    } while (nmodif);
  
    //MAJ des refs
    fit = M_faceIter(mesh);
    while ((pface = FIter_next(fit))) {
      assert(!pface->deleted);
      int num = F_numRegions(pface);
      if(!num) {
        printf("face %p dim %d\n",(void*)pface,GEN_type(pface->g));
        M_removeFace(mesh,pface); 
        continue;
      }
      int dim = GEN_type(pface->g); 
      if(dim == 2) {
        int num = F_numRegions(pface);
        if(num == 2) {
          //puts("enlever la ref");
          pRegion pr = F_region(pface,0);
          //	pGEntity pg  = EN_whatIn(pface);
          pGEntity pgr = EN_whatIn(pr);
          int tagr = GEN_tag(pgr);
          int dimr = GEN_type(pgr); 
          if(dimr==3) {
            pface->g = (pGEntity) GM_regionByTag(mesh->model,tagr);   
            assert(   GEN_type(pface->g)==3);
          } else {
            printf("----pbs faces**** %d\n",dim);
          }

        } else if(!num) {
          M_removeFace(mesh,pface); 
        }
        if(!(F_numRegions(pface)==1 || GEN_type(pface->g)!=2)) {
          printf("face %p %d %d %d\n",(void*)pface,F_numRegions(pface),GEN_type(pface->g),pface->deleted);
        }
        assert(F_numRegions(pface)==1 || GEN_type(pface->g)!=2 || pface->deleted);
      }
    }
    FIter_delete(fit);
    EIter eit = M_edgeIter(mesh);
    pEdge ped;
    while ((ped = EIter_next(eit))) {
      int num = E_numFaces(ped);
      if(!num) {
        printf("edge dim %d\n",GEN_type(ped->g));
        M_removeEdge(mesh,ped); 
        continue;
      }

    }
    EIter_delete(eit);
  }


  void GroupPeriodicTetra(pMesh mesh, MDB_DataExchanger &de,
                          MDB_DataExchangerPeriodic &deperiodic) {

    pMeshDataId tagMove = MD_newMeshDataId("TagMovePeriodic");
    pMeshDataId tagElt  = MD_newMeshDataId("EltDestination"); //dest = (int - 1)

    int maxiter = 5;
    int iter = 1;                                       
    while(MarkGroupPeriodicTets(mesh,tagElt,tagMove) && iter++ < maxiter) { 
      printf("group iter %d\n",iter-1);
      /* 1) marked vertex*/
      MarkPeriodicEltVertex(mesh,tagElt,tagMove);

      /* 2) move elements*/                                         
      MovePeriodicTetrasInverse(mesh,tagElt,tagMove,de,deperiodic);
  
      /* 3) delete vertex*/
      VertexDelete(mesh,tagMove,tagMove); 
  
    } 
    MD_deleteMeshDataId(tagMove);
    MD_deleteMeshDataId(tagElt);  
  
    return;
  }

} // End of namespace MAd
