// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

// GC: this should be moved to Mesh/

#include "MeshDataBase.h"
#include "MeshDataBaseParallelInterface.h"

#include "Constraint.h"
#include "ModelConstraintManager.h"

static unsigned int CONSTRAIN_MARK_ID = 1687354269;

namespace MAd {

  // -------------------------------------------------------------------
  void EN_constrain(pEntity pE) 
  {
    pE->attachInt(CONSTRAIN_MARK_ID,1);
  }

  // -------------------------------------------------------------------
  void EN_unconstrain(pEntity pE)
  {
    pE->attachInt(CONSTRAIN_MARK_ID,0);
  }
 
  // -------------------------------------------------------------------
  bool EN_constrained(pEntity pE)
  {
    if ( ModelConstraintManagerSgl::instance().constrained(EN_whatIn(pE)) ) {
      return true;
    }

    return (bool) pE->getAttachedInt(CONSTRAIN_MARK_ID);
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void DeleteConstraint(pMesh mesh) {
    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      EN_unconstrain((pEntity) pv);  
    }
    VIter_delete(vit);	
  
    EIter eit = M_edgeIter(mesh);
    pEdge ped;  
    while ((ped = EIter_next(eit))) {
      EN_unconstrain((pEntity) ped);  
    }
    EIter_delete(eit);	

    FIter fit = M_faceIter(mesh);
    pFace pface;  
    while ((pface = FIter_next(fit))) {
      EN_unconstrain((pEntity) pface);  
    }
    FIter_delete(fit);	

    RIter rit = M_regionIter(mesh);
    pRegion pregion;  
    while ((pregion = RIter_next(rit))) {
      EN_unconstrain((pEntity) pregion);  
    }
    RIter_delete(rit);	
  }

  // -------------------------------------------------------------------
#ifdef PARALLEL
  void DeleteParallelConstraint(pMesh mesh) {

    VIter vit = M_vertexIter(mesh);
    while ( pVertex pv = VIter_next(vit) ) {
      if ( V_isInterface(pv) ) EN_unconstrain((pEntity) pv);
    }
    VIter_delete(vit);

    EIter eit = M_edgeIter(mesh);
    while ( pEdge pe = EIter_next(eit) ) {
      // if ( E_isInterface(mesh,pe) ) {
      if ( E_isInterface(pe) ) {
        EN_unconstrain((pEntity) pe);
      }
    }
    EIter_delete(eit);
    int dim = (mesh->tets.empty()) ? 2 : 3;
  
    if ( dim == 3 ) {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        // if ( F_isInterface(mesh,pf) ) {
        if ( F_isInterface(pf) ) {
          EN_unconstrain((pEntity) pf);
        }
      }
      FIter_delete(fit);
    }
  }
  void UpdateParallelConstraint(pMesh mesh) {

    DeleteParallelConstraint(mesh);
    //   DeleteConstraint(mesh);
    int dim = (mesh->tets.empty()) ? 2 : 3;
  
    if ( dim == 3 ) {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        // if ( F_isInterface(mesh,pf) ) {
        if ( F_isInterface(pf) ) {
          EN_constrain((pEntity) pf);
        }
      }
      FIter_delete(fit);
    }
    EIter eit = M_edgeIter(mesh);
    while ( pEdge pe = EIter_next(eit) ) {
      // if ( E_isInterface(mesh,pe) ) {
      if ( E_isInterface(pe) ) {
        EN_constrain((pEntity) pe);
      }
    }
    EIter_delete(eit);
    VIter vit = M_vertexIter(mesh);
    while ( pVertex pv = VIter_next(vit) ) {
      if ( V_isInterface(pv) ) EN_constrain((pEntity) pv);
    }
    VIter_delete(vit);
  }
#endif

  // -------------------------------------------------------------------
  void UpdatePeriodicConstraint2d(pMesh mesh)
  {
    pMeshDataId tagPeriod  = MD_lookupMeshDataId("PeriodicPoint");

    FIter fit = M_faceIter(mesh);
    pFace pface;  
    while ((pface = FIter_next(fit))) {
      int ncount = 0;
      pEdge pe1 = F_edge(pface,0);
      pVertex p1 = E_vertex(pe1,0);
      pVertex p2 = E_vertex(pe1,1);
      void *temp_ptr1,*temp_ptr2;
      int isPeriodic1  = EN_getDataPtr((pEntity) p1 , tagPeriod, &temp_ptr1);
      int isPeriodic2  = EN_getDataPtr((pEntity) p2 , tagPeriod, &temp_ptr2);
      if(isPeriodic1) EN_constrain((pEntity) p1); 
      if(isPeriodic2) EN_constrain((pEntity) p2);
      if(isPeriodic1 && isPeriodic2){ 
        /*est-ce que l'arete image de celle-ci existe ?*/
        /*ie est-ce l'arete Img(P1)-Img(P2) existe et est differente de ped ?*/
        const std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;
        const std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
        int nump = EN_id(p2); 
        unsigned int k=0;
        int exist=-1;
        for(k=0; k<(*recup).size() ; k++) { 
          /*edge p1p2 = p1 Img(p1) ?*/
          pVertex n1 = (*recup)[k].second;
          if ( nump == EN_id(n1) ) {
            continue;         		
          } 
          /*edge = Img(p1)Img(p2) existe ?*/  
          unsigned int kk=0;
          std::vector<int> transfo1 = (*recup)[k].first;
          for(kk=0; kk<(*recup2).size() ; kk++) { 
            std::vector<int> transfo2 = (*recup2)[kk].first;
            assert(transfo1.size()==transfo2.size());
            unsigned int kt = 0;
            for(kt = 0 ; kt < transfo1.size(); kt++) {
              if(transfo1[kt] != transfo2[kt]) break;
            }
            if(kt!=transfo1.size()) continue;
            pVertex n2 = (*recup2)[kk].second;   
            if ( !E_exist(n1,n2) && !E_exist(n2,n1) ) { 
              exist = 0;         		
            } else {
              exist = 1;
            } 
          }
        }
        if(exist==1) {      
          ncount++;
          EN_constrain((pEntity) pe1);
        }
      }
      pEdge pe2 = F_edge(pface,1);
      p1 = E_vertex(pe2,0);
      p2 = E_vertex(pe2,1);
      isPeriodic1  = EN_getDataPtr((pEntity) p1 , tagPeriod, &temp_ptr1);
      isPeriodic2  = EN_getDataPtr((pEntity) p2 , tagPeriod, &temp_ptr2);
      if(isPeriodic1) EN_constrain((pEntity) p1); 
      if(isPeriodic2) EN_constrain((pEntity) p2);
      if(isPeriodic1 && isPeriodic2){ 
        /*est-ce que l'arete image de celle-ci existe ?*/
        /*ie est-ce l'arete Img(P1)-Img(P2) existe et est differente de ped ?*/
        const std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;
        const std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
        int nump = EN_id(p2); 
        unsigned int k=0;
        int exist=-1;
        for(k=0; k<(*recup).size() ; k++) {
          /*edge p1p2 = p1 Img(p1) ?*/
          pVertex n1 = (*recup)[k].second;
          if ( nump == EN_id(n1) ){
            continue;         		
          } 
          /*edge = Img(p1)Img(p2) existe ?*/  
          unsigned int kk=0;
          std::vector<int> transfo1 = (*recup)[k].first;
          for(kk=0; kk<(*recup2).size() ; kk++) { 
            std::vector<int> transfo2 = (*recup2)[kk].first;
            assert(transfo1.size()==transfo2.size());
            unsigned int kt = 0;
            for(kt = 0 ; kt < transfo1.size(); kt++) {
              if(transfo1[kt] != transfo2[kt]) break;
            }
            if(kt!=transfo1.size()) continue;
            pVertex n2 = (*recup2)[kk].second;   
            if ( !E_exist(n1,n2) && !E_exist(n2,n1) ) { 
              exist = 0;         		
            } else {
              exist = 1;
            } 
          }
        }
        if(exist==1) {      
          ncount++;
          EN_constrain((pEntity) pe2); 	
        }
      }	
      pEdge pe3 = F_edge(pface,2);
      p1 = E_vertex(pe3,0);
      p2 = E_vertex(pe3,1);
      isPeriodic1  = EN_getDataPtr((pEntity) p1 , tagPeriod, &temp_ptr1);
      isPeriodic2  = EN_getDataPtr((pEntity) p2 , tagPeriod, &temp_ptr2);
      if(isPeriodic1) EN_constrain((pEntity) p1); 
      if(isPeriodic2) EN_constrain((pEntity) p2);
      if(isPeriodic1 && isPeriodic2){ 
        /*est-ce que l'arete image de celle-ci existe ?*/
        /*ie est-ce l'arete Img(P1)-Img(P2) existe et est differente de ped ?*/
        const std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;
        const std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
        int nump = EN_id(p2); 
        unsigned int k=0,exist=0;
        for(k=0; k<(*recup).size() ; k++) {
          /*edge p1p2 = p1 Img(p1) ?*/
          pVertex n1 = (*recup)[k].second;
          if ( nump == EN_id(n1) ){
            continue;         		
          } 
          /*edge = Img(p1)Img(p2) existe ?*/  
          unsigned int kk=0;
          std::vector<int> transfo1 = (*recup)[k].first;
          for(kk=0; kk<(*recup2).size() ; kk++) { 
            std::vector<int> transfo2 = (*recup2)[kk].first;
            assert(transfo1.size()==transfo2.size());
            unsigned int kt = 0;
            for(kt = 0 ; kt < transfo1.size(); kt++) {
              if(transfo1[kt] != transfo2[kt]) break;
            }
            if(kt!=transfo1.size()) continue;
            pVertex n2 = (*recup2)[kk].second;   
            if ( E_exist(n1,n2) || E_exist(n2,n1) ) { 
              exist = 1;
            } 
          }
        }
        if(exist==1) {      
          ncount++;
          EN_constrain((pEntity) pe3);
        }
      }
      if(ncount==3)  EN_constrain((pEntity) pface);
  
    }
    FIter_delete(fit);
  }

  // -------------------------------------------------------------------
  void UpdatePeriodicConstraint3d(pMesh mesh)
  {
    pMeshDataId tagPeriod  = MD_lookupMeshDataId("PeriodicPoint");

    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      int ncountr = 0;
      pFace pface[4];
      for(int i=0 ; i<4 ; i++) { 
        int nbdry = 1;
        pface[i] = R_face(pr,i);  
        /*Img(pface[i]) exist ?*/
        pVertex nod[3];
        int isPeriodic[3];
        void *tmp[3];
        unsigned int j=0;
        for(j=0 ;  j<3 ; j++) {
          nod[j] = F_vertex(pface[i],j);  
          isPeriodic[j]  = EN_getDataPtr((pEntity) nod[j] , tagPeriod, &tmp[j]);
          if(!isPeriodic[j]) {
            nbdry=0;
            continue;
          }
          EN_constrain((pEntity) nod[j]); 
        }             
        if(nbdry) {
          /*find Img(nod[j])*/ 
          const std::vector<std::pair<std::vector<int> , pVertex> > *recup  = (std::vector<std::pair<std::vector<int> , pVertex> > *) tmp[0];
          const std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) tmp[1];
          const std::vector<std::pair<std::vector<int> , pVertex> > *recup3 = (std::vector<std::pair<std::vector<int> , pVertex> > *) tmp[2];
          unsigned int k;
          int exist=-1,nfound=0;
          for(k=0; k<(*recup).size() ; k++) {
            pVertex n1 = (*recup)[k].second;
            /*ATTENTION CA PEUT ARRIVER DANS LES CAS 3_PERIODIQUES!!!if ( EN_id((pEntity) nod[1]) == EN_id((pEntity) n1) || EN_id((pEntity) nod[2]) == EN_id((pEntity) n1) ){  
              continue;         		
              } */
            unsigned int kk=0;
            std::vector<int> transfo1 = (*recup)[k].first;
            for(kk=0; kk<(*recup2).size() ; kk++) { 
              std::vector<int> transfo2 = (*recup2)[kk].first;
              assert(transfo1.size()==transfo2.size());
              unsigned int kt = 0;
              for(kt = 0 ; kt < transfo1.size(); kt++) {
                if(transfo1[kt] != transfo2[kt]) break;
              }
              if(kt!=transfo1.size()) continue; 
              if(!nfound) nfound = 1;
              pVertex n2 = (*recup2)[kk].second;   
              unsigned int kkk=0;
              for(kkk=0; kkk<(*recup3).size() ; kkk++) { 
                std::vector<int> transfo2 = (*recup3)[kkk].first;
                unsigned int kt = 0;
                for(kt = 0 ; kt < transfo1.size(); kt++) {
                  if(transfo1[kt] != transfo2[kt]) break;
                }
                if(kt!=transfo1.size()) continue;
                pVertex n3 = (*recup3)[kkk].second;   
                nfound = 2;
                if ( F_exist(n1,n2,n3,0) || F_exist(n1,n3,n2,0) ||
                     F_exist(n2,n1,n3,0) || F_exist(n2,n3,n1,0) || 
                     F_exist(n3,n2,n1,0) || F_exist(n3,n1,n2,0) ) { 
                  exist = 1;
                } 
              }/*end kkk*/    
            }/*end kk*/
          }/*end k*/ 
          if(!nfound) nbdry=0;//printf("on n'a pas trouve la face\n"); 
          if(nfound==1) nbdry = 1;//printf("3 ou 2 ? %d\n",nfound);
          if(nfound==2) nbdry = 0;//printf("3 ou 2 ? %d\n",nfound);
	
          if(exist==1) {      
            ncountr++;
            EN_constrain((pEntity) pface[i]);
            for(int j=0 ;  j<3 ; j++) {
              EN_constrain((pEntity) F_edge(pface[i],j));  
              /* pVertex p1 = F_edge(pface[i],j)->p1;
                 pVertex p2 = F_edge(pface[i],j)->p2;
                 double len = (p1->X-p2->X)*(p1->X-p2->X) + (p1->Y-p2->Y)*(p1->Y-p2->Y) +(p1->Z-p2->Z)*(p1->Z-p2->Z);
                 printf("edge %d %d : %e\n",p1->iD,p2->iD,len); */
            }
            continue;
          } else{
            /* if(F_numRegions(pface[i]) != 2){ printf("pbs constraint %d %d %d\n",EN_id((pEntity) nod[0])
               ,EN_id((pEntity) nod[1]),EN_id((pEntity) nod[2]));exit(0);  }   */
          }
        }/*end nbdry*/                    
        /*face not constraint*/
        int ncount = 0;      
        for(int j=0 ; j<3 ; j++) {
          pEdge ped = F_edge(pface[i],j);
          pVertex p1 = E_vertex(ped,0);
          pVertex p2 = E_vertex(ped,1);
          void *temp_ptr1,*temp_ptr2;
          int isPeriodic1  = EN_getDataPtr((pEntity) p1 , tagPeriod, &temp_ptr1);
          int isPeriodic2  = EN_getDataPtr((pEntity) p2 , tagPeriod, &temp_ptr2);
          if((isPeriodic1 && isPeriodic2)) {
            const std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr1;
            const std::vector<std::pair<std::vector<int> , pVertex> > *recup2 = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr2;
            //int nump = EN_id(p2);
            //int nfound = 0;
            unsigned int k=0;
            int exist=-1;
            for(k=0; k<(*recup).size() ; k++) {
              /*edge p1p2 = p1 Img(p1) ?*/
              pVertex n1 = (*recup)[k].second;  
              /*if ( nump == EN_id(n1) ){
                continue;         		
                }*/ 
              /*edge = Img(p1)Img(p2) existe ?*/  
              unsigned int kk=0;
              std::vector<int> transfo1 = (*recup)[k].first;
              for(kk=0; kk<(*recup2).size() ; kk++) { 
                std::vector<int> transfo2 = (*recup2)[kk].first;
                assert(transfo1.size()==transfo2.size());
                unsigned int kt = 0;
                for(kt = 0 ; kt < transfo1.size(); kt++) {
                  if(transfo1[kt] != transfo2[kt]) break;
                }
                if(kt!=transfo1.size()) continue;
                //nfound = 1;
                pVertex n2 = (*recup2)[kk].second;   
                if ( E_exist(n1,n2) || E_exist(n2,n1) ) { 
                  exist = 1;
                } 
              }
            }/*end k*/
		   
            if(exist==1) {
              ncount++; 
              EN_constrain((pEntity) ped);
            }  
          }
        }
      }
      if(ncountr==4) EN_constrain((pEntity) pr);
    }
    RIter_delete(rit); 
  
  }

}

//#endif
