// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Gaetan Compere, Olivier Pierard, Jean-Francois Remacle
// -------------------------------------------------------------------

#include "OperatorTools.h"
#include "CallbackManager.h"
#include "MAdMessage.h"
#include "Constraint.h"

#include <map>
#include <assert.h>

namespace MAd {

  // -------------------------------------------------------------------
  bool V_setPosition(pVertex vertex, double target[3])
  {
    if( EN_constrained((pEntity)vertex) ) return false;
    CallBackManagerSgl::instance().callCallBackMoves(vertex,target);
    pPoint point = V_point(vertex);
    P_setPos(point,target[0],target[1],target[2]);
    return true;
  }

  // -------------------------------------------------------------------
  void E_collapse(pMesh mesh, 
                  pEdge edgeDel, 
                  pVertex vDel, 
                  pVertex vTgt)
  {
    pPList eRegions = E_regions(edgeDel);

    // --- 2D collapse case ---
    if( PList_size(eRegions) == 0 ) {
      PList_delete(eRegions);
      E_collapseOnGFace(mesh,edgeDel,vDel,vTgt);
      return;
    }

    // --- build the new mesh ---

    pVertex verts[4];
    pPList newRegions = PList_new();
    pPList vDelRegions = V_regions(vDel);
    void * temp = NULL;
    while ( pRegion region = (pRegion)PList_next(vDelRegions,&temp) )
      {
        if( PList_inList(eRegions,(pEntity)region) ) continue;
        
        pFace fExt = R_vtOpFc(region,vDel);
        
        // create the region
        pGEntity rGEntity = (pGEntity)R_whatIn(region);
        verts[0] = vTgt;
        for (int iV = 0; iV<3; iV++) verts[iV+1] = F_vertex(fExt,iV);
        pRegion newR = M_createR (mesh, 4, verts, rGEntity);
        PList_append(newRegions, (pEntity)newR);
      }
    PList_delete(eRegions);

    // --- classify new boundary entities ---
    int gDim = E_whatInType(edgeDel);

    // Classify edges on lines
    if ( gDim == 1 )
      {
        pPList vDelEdges = V_edges(vDel);
        temp = NULL;
        pEdge pe;
        while ( ( pe = (pEdge)PList_next(vDelEdges,&temp) ) ) {
          if ( pe == edgeDel ) continue;
          if ( E_whatInType(pe) == 1 ) {
            pGEntity gEdge2 = E_whatIn(pe);
            pVertex oppV = E_otherVertex(pe,vDel);
            pEdge newLineEdge = E_exist(vTgt,oppV);
            E_setWhatIn(newLineEdge,gEdge2);
          }
        }
        PList_delete(vDelEdges);
      }

    // classify faces on surfaces
    if ( gDim <= 2 )
      {
        pGEntity gFace;
        pEdge oppE;
        pFace newSurfaceFace;

        pPList vDelFaces = V_faces(vDel);
        temp = NULL;
        pFace pf;
        while ( ( pf = (pFace)PList_next(vDelFaces,&temp) ) ) {
          if ( F_inClosure(pf,(pEntity)edgeDel) ) continue;
          if ( F_whatInType(pf) == 2 ) {
            gFace          = F_whatIn(pf);
            oppE           = F_vtOpEd(pf,vDel);
            newSurfaceFace = F_exist(oppE,vTgt);
            F_setWhatIn(newSurfaceFace,gFace);
          }
        }
        PList_delete(vDelFaces);
      }

    // --- call callback functions ---
    CallBackManagerSgl::instance().callCallBacks(vDelRegions,
                                                 newRegions,
                                                 MAd_ECOLLAPSE,
                                                 (pEntity)vDel);
    PList_delete(newRegions);

    // --- delete old cavity ---
    temp = NULL;
    while ( pRegion region = (pRegion)PList_next(vDelRegions,&temp) ) {
      M_removeRegion(mesh,region);
    }
    PList_delete(vDelRegions);

    pPList vFaces = V_faces(vDel);
    temp = NULL;
    while ( pFace face = (pFace)PList_next(vFaces,&temp) ) {
      M_removeFace(mesh,face);
    }
    PList_delete(vFaces);

    pPList vEdges = V_edges(vDel);
    temp = NULL;
    while ( pEdge edge = (pEdge)PList_next(vEdges,&temp) ) {
      M_removeEdge(mesh,edge);
    }
    PList_delete(vEdges);

    M_removeVertex(mesh,vDel);

  }
    

    
    /*
// This is the old way to make a collapse

    pPList eRegions = E_regions(edgeDel);

    // --- 2D collapse case ---
    if( PList_size(eRegions) == 0 ) {
      PList_delete(eRegions);
      E_collapseOnGFace(mesh,edgeDel,vDel,vTgt);
      return;
    }

    // --- list faces to be merged ---

    std::map<pFace,pFace> fOldToNew; // deleted face -> corresponding new (or target) face
    std::map<pFace,bool>  fDirs;     // if direction of (deleted) face has to be changed

    void * temp1 = 0;
    while ( pRegion region = (pRegion)PList_next(eRegions,&temp1) ) {
      
      pFace fDel, fTgt;
      
      pPList rFaces = R_faces(region);
      void * temp2 = 0;
      while ( pFace face = (pFace)PList_next(rFaces, &temp2) ) {
        if (!F_inClosure(face,(pEntity)edgeDel)) {
          if (F_inClosure(face,(pEntity)vDel)) fDel = face;  // deleted face
          else                                 fTgt = face;  // target face
        }
      }
      PList_delete(rFaces);

      fOldToNew[fDel] = fTgt;
      
      // Choose the direction and classification of target face
      pEdge edge = F_vtOpEd(fDel,vDel);
      if( F_whatInType(fDel) < F_whatInType(fTgt) ) {
        F_setWhatIn( fTgt, F_whatIn(fDel) );
        if( F_dirUsingEdge(fTgt,edge) != F_dirUsingEdge(fDel,edge) ) {
          F_chDir(fTgt);
        }
      }
      else {
        if( F_dirUsingEdge(fTgt,edge) == F_dirUsingEdge(fDel,edge) ) {
          fDirs[fDel] = true;
        } 
        else {
          fDirs[fDel] = false;
        }
      }
    }

    // --- list edges to be merged ---

    std::map<pEdge,pEdge> eOldToNew; // deleted edge -> corresponding new (or target) edge
    std::map<pEdge,bool>  eDirs;     // if direction of (deleted) edge has to be changed

    pPList eFaces = E_faces(edgeDel);
    void * temp3 = 0;
    while ( pFace face = (pFace)PList_next(eFaces,&temp3) ) {
      
      pEdge eDel, eTgt;

      pPList fEdges = F_edges(face);
      void * temp4 = 0;
      while ( pEdge edge = (pEdge)PList_next(fEdges, &temp4) ) {
        if ( edge != edgeDel ) {
          if (E_inClosure(edge,(pEntity)vDel)) eDel = edge;  // deleted edge
          else                                 eTgt = edge;  // target edge
        }
      }
      PList_delete(fEdges);

      eOldToNew[eDel] = eTgt;
      
      // Choose the classification of target edge
      if( E_whatInType(eDel) < E_whatInType(eTgt) ) {
        E_setWhatIn( eTgt, E_whatIn(eDel) );
      }

      // Choose the direction of target edge
      if( E_vertex(eDel,0) == E_vertex(eTgt,0) ||
          E_vertex(eDel,1) == E_vertex(eTgt,1) ) {
        eDirs[eDel] = true;
      }
      else {
        eDirs[eDel] = false;
      }
    }
    PList_delete(eFaces);

    // --- build the new mesh ---

    pPList newRegions = PList_new();
    pPList vDelRegions = V_regions(vDel);
    void * temp5 = 0;
    while ( pRegion region = (pRegion)PList_next(vDelRegions,&temp5) ) {

      if( PList_inList(eRegions,(pEntity)region) ) continue;

      // need to create a new region
      pFace rFaces[4];

      pFace fExt = R_vtOpFc(region,vDel);

      for(int iF=0; iF<4; iF++) {
        
        pFace face = R_face(region,iF);

        if( face == fExt ) {
          rFaces[iF] = face;
          continue;
        }
        
        if ( fOldToNew.find(face) != fOldToNew.end() ) {
          rFaces[iF] = fOldToNew.find(face)->second;
        }
        else {
          // need to create a new face
          pEdge fEdges[3];

          for(int iE=0; iE<3; iE++ ) {

            pEdge edge = F_edge(face,iE);

            if( F_inClosure(fExt,(pEntity)edge) ) {
              fEdges[iE] = edge;
              continue;
            }
        
            if ( eOldToNew.find(edge) != eOldToNew.end() ) {
              fEdges[iE] = eOldToNew.find(edge)->second;
            }
            else {
              // need to create a new edge

              pVertex eVertices[2];
              for(int iV=0; iV<2; iV++) {
                pVertex pV = E_vertex(edge,iV);
                if( pV == vDel ) eVertices[iV] = vTgt;
                else             eVertices[iV] = pV;
              }
          
              // create the edge
              pGEntity eGEntity = E_whatIn(edge);
              fEdges[iE] = M_createE(mesh, eVertices[0], eVertices[1], eGEntity);
              eOldToNew[edge] = fEdges[iE];
            }
          }
      
          // create the face
          pGEntity fGEntity = F_whatIn(face);
          rFaces[iF] = M_createF(mesh, 3, fEdges, fGEntity);
          fOldToNew[face] = rFaces[iF];
        }
      }

      // create the region
      pGEntity rGEntity = (pGEntity)R_whatIn(region);
      pRegion newRegion = M_createR(mesh, 4, rFaces, rGEntity);
      PList_append(newRegions, (pEntity)newRegion);
    }

    PList_delete(eRegions);

    // --- call callback functions ---
    CallBackManagerSgl::instance().callCallBacks(vDelRegions,
                                                 newRegions,
                                                 MAd_ECOLLAPSE,
                                                 (pEntity)vDel);
    PList_delete(newRegions);

    // --- delete old cavity ---
    void * temp = 0;
    while ( pRegion region = (pRegion)PList_next(vDelRegions,&temp) ) {
      M_removeRegion(mesh,region);
    }
    PList_delete(vDelRegions);

    pPList vFaces = V_faces(vDel);
    temp = 0;
    while ( pFace face = (pFace)PList_next(vFaces,&temp) ) {
      M_removeFace(mesh,face);
    }
    PList_delete(vFaces);

    pPList vEdges = V_edges(vDel);
    temp = 0;
    while ( pEdge edge = (pEdge)PList_next(vEdges,&temp) ) {
      M_removeEdge(mesh,edge);
    }
    PList_delete(vEdges);

    M_removeVertex(mesh,vDel);

    */

  // -------------------------------------------------------------------
  void E_collapseOnGFace(pMesh mesh, 
                         pEdge edgeDel, 
                         pVertex vDel, 
                         pVertex vTgt)
  {
    // --- list edges to be merged ---

    std::map<pEdge,pEdge> eOldToNew; // deleted edge -> corresponding new (or target) edge
    std::map<pEdge,bool>  eDirs;     // if direction of (deleted) edge has to be changed

    pPList eFaces = E_faces(edgeDel);
    void * temp1 = 0;
    while ( pFace face = (pFace)PList_next(eFaces,&temp1) ) {
      
      pEdge eDel, eTgt;

      pPList fEdges = F_edges(face);
      void * temp2 = 0;
      while ( pEdge edge = (pEdge)PList_next(fEdges, &temp2) ) {
        if ( edge != edgeDel ) {
          if (E_inClosure(edge,(pEntity)vDel)) eDel = edge;  // deleted edge
          else                                 eTgt = edge;  // target edge
        }
      }
      PList_delete(fEdges);

      eOldToNew[eDel] = eTgt;
      
      // Choose the classification of target edge
      if( E_whatInType(eDel) < E_whatInType(eTgt) ) {
        E_setWhatIn( eTgt, E_whatIn(eDel) );
      }

      // Choose the direction of target edge
      if( E_vertex(eDel,0) == E_vertex(eTgt,0) ||
          E_vertex(eDel,1) == E_vertex(eTgt,1) ) {
        eDirs[eDel] = true;
      }
      else {
        eDirs[eDel] = false;
      }
    }

    // --- build the new mesh ---

    pPList newFaces = PList_new();
    pPList vDelFaces = V_faces(vDel);
    void * temp3 = 0;
    while ( pFace face = (pFace)PList_next(vDelFaces,&temp3) ) {

      if( PList_inList(eFaces,(pEntity)face) ) continue;

      // need to create a new face
      pEdge fEdges[3];

      pEdge eExt = F_vtOpEd(face,vDel);

      for(int iE=0; iE<3; iE++ ) {

        pEdge edge = F_edge(face,iE);

        if( edge == eExt ) {
          fEdges[iE] = edge;
          continue;
        }
        
        if ( eOldToNew.find(edge) != eOldToNew.end() ) {
          fEdges[iE] = eOldToNew.find(edge)->second;
        }
        else {
          // need to create a new edge

          pVertex eVertices[2];
          for(int iV=0; iV<2; iV++) {
            pVertex pV = E_vertex(edge,iV);
            if( pV == vDel ) eVertices[iV] = vTgt;
            else             eVertices[iV] = pV;
          }
          
          // create the edge
          pGEntity eGEntity = E_whatIn(edge);
          fEdges[iE] = M_createE(mesh, eVertices[0], eVertices[1], eGEntity);
          eOldToNew[edge] = fEdges[iE];
        }
      }
      
      // create the face
      pGEntity fGEntity = F_whatIn(face);
      pFace newFace = M_createF(mesh, 3, fEdges, fGEntity);
      PList_append(newFaces,(pEntity)newFace);
    }

    PList_delete(eFaces);

    // --- call callback functions ---
    CallBackManagerSgl::instance().callCallBacks(vDelFaces,
                                                 newFaces,
                                                 MAd_ECOLLAPSE,
                                                 (pEntity)vDel);
    PList_delete(newFaces);

    // --- delete old cavity ---
    void * temp = 0;
    while ( pFace face = (pFace)PList_next(vDelFaces,&temp) ) {
      M_removeFace(mesh,face);
    }
    PList_delete(vDelFaces);

    pPList vEdges = V_edges(vDel);
    temp = 0;
    while ( pEdge edge = (pEdge)PList_next(vEdges,&temp) ) {
      M_removeEdge(mesh,edge);
    }
    PList_delete(vEdges);

    M_removeVertex(mesh,vDel);
  }


  // -------------------------------------------------------------------
  pVertex E_split(pMesh mesh, pEdge edge, double xyz[3], double t)
  {
    // --- build the vertex ---

    pVertex newV = NULL; // the new vertex

    pGEntity edgeGE = E_whatIn(edge); // its classification

    double u[2][2];
    u[0][0] = -2.; u[0][1] = -2.; u[1][0] = -2.; u[1][1] = -2.;
    if ( E_params(edge,u) ) {
      double uV[2];
      if ( GEN_type(edgeGE) == 2 ) {
        GF_centerOnGeodesic( (pGFace)edgeGE, t, u, uV );
      }
      else {
        for (int iP=0; iP<2; iP++)  uV[iP] = (1.-t) * u[0][iP] + t * u[1][iP];
      }
      newV = M_createVP(mesh, xyz[0], xyz[1], xyz[2], uV[0], uV[1], -1, edgeGE);
    }
    else {
      newV = M_createV(mesh, xyz[0], xyz[1], xyz[2], -1, edgeGE);
    }

    // --- build the edges ---
    pPList newEdges = PList_new();
    pEdge e1 = M_createE(mesh, E_vertex(edge,0), newV, edgeGE);
    PList_append(newEdges,(pEntity)e1);
    pEdge e2 = M_createE(mesh, newV, E_vertex(edge,1), edgeGE);
    PList_append(newEdges,(pEntity)e2);
    
    // --- build the faces bordered by new edges ---
    pPList newFaces = PList_new();
    pPList eFaces = E_faces(edge);
    void * temp = NULL;
    while( pFace face = (pFace)PList_next(eFaces,&temp) ) {

      pGEntity fGE = F_whatIn(face);
      pEdge newE = M_createE(mesh, newV, F_edOpVt(face,edge), fGE);
      
      int iEdge = 0;
      for( iEdge=0; iEdge< F_numEdges(face); iEdge++) {
        if( F_edge(face,iEdge) == edge ) break ;
      }
      pEdge fEdges[3];
      fEdges[0] = edge ;
      fEdges[1] = F_edge(face,(iEdge+1)%3);
      fEdges[2] = F_edge(face,(iEdge+2)%3);

      pFace newF1, newF2;
      if( F_edgeDir(face,iEdge) ) {

        pEdge newFEdges1[3] = {e1, newE, fEdges[2]};
        newF1 = M_createF(mesh, 3, newFEdges1, fGE);
        EN_attachDataP((pEntity)fEdges[2],"_EdgeSplitMarker_",(void *)newF1);

        pEdge newFEdges2[3] = {e2, fEdges[1], newE};
        newF2 = M_createF(mesh, 3, newFEdges2, fGE);
        EN_attachDataP((pEntity)fEdges[1],"_EdgeSplitMarker_",(void *)newF2);
      }
      else {
        pEdge newFEdges1[3] = {e1, fEdges[1], newE};
        newF1 = M_createF(mesh, 3, newFEdges1, fGE);
        EN_attachDataP((pEntity)fEdges[1],"_EdgeSplitMarker_",(void *)newF1); 
        
        pEdge newFEdges2[3] = {e2, newE, fEdges[2]};
        newF2 = M_createF(mesh, 3, newFEdges2, fGE);
        EN_attachDataP((pEntity)fEdges[2],"_EdgeSplitMarker_",(void *)newF2);
      }
      EN_attachDataP((pEntity)face,"_EdgeSplitMarker_",(void *)newE);
      EN_attachDataP((pEntity)newF1,"_EdgeSplitMarker_",face);
      EN_attachDataP((pEntity)newF2,"_EdgeSplitMarker_",face);
      
      PList_append(newFaces,(pEntity)newF1);
      PList_append(newFaces,(pEntity)newF2);
    }
    
    pPList delRegs = PList_new();
    pPList eRegs = E_regions(edge);
    if (eRegs) {

      temp = NULL;
      while( pRegion region = (pRegion)PList_next(eRegs, &temp) ) {

        pGEntity rGE = (pGEntity)R_whatIn(region);

        pFace extFaces[2];

        // --- build the new face between the two new regions ---
        pEdge fEdges[3];
        for(int iF=0, eCount=0, fCount=0; iF<R_numFaces(region); iF++) {
          pFace face = R_face(region,iF);
          if( F_inClosure(face,(pEntity)edge) ) {
            fEdges[eCount++] = (pEdge)EN_dataP((pEntity)face,"_EdgeSplitMarker_");
          }
          else {
            extFaces[fCount++] = face;
          }
        }
        pEdge oppE = R_gtOppEdg(region,edge);
        fEdges[2] = oppE;
        //         if( E_vertex(fEdges[0],1) != E_vertex(fEdges[1],0) && 
        //             E_vertex(fEdges[0],1) != E_vertex(fEdges[1],1) )  {
        //           pEdge tmpE = fEdges[1] ;
        //           fEdges[1] = fEdges[2];
        //           fEdges[2] = tmpE;
        //         }
        pFace newF = M_createF(mesh, 3, fEdges, rGE);

        // -- build the regions ---
        pFace rFaces[4];
        
        rFaces[0] = extFaces[0] ;
        rFaces[1] = newF;
        int fIndex = 2;
        for(int i=0; i<F_numEdges(extFaces[0]); i++) {
          pEdge tmpE = F_edge(extFaces[0],i);
          if( tmpE != oppE ) {
            rFaces[fIndex] = (pFace)EN_dataP((pEntity)tmpE,"_EdgeSplitMarker_");
            fIndex++;
          }
        }
        M_createR(mesh, 4, rFaces, rGE);

        rFaces[0] = extFaces[1] ;
        rFaces[1] = newF;
        fIndex = 2;
        for( int i=0; i < F_numEdges(extFaces[1]); i++ ) {
          pEdge tmpE = F_edge(extFaces[1],i);
          if( tmpE != oppE ) {
            rFaces[fIndex] = (pFace)EN_dataP((pEntity)tmpE,"_EdgeSplitMarker_");
            fIndex++;
          }
        }
        M_createR(mesh, 4, rFaces, rGE);

        PList_append(delRegs,(pEntity)region);
      }
    }
    PList_delete(eRegs);
    
    pPList delEdge = PList_new();
    PList_append(delEdge,(pEntity)edge);
    CallBackManagerSgl::instance().callCallBacks(delEdge, newEdges, 
                                                 MAd_ESPLIT, (pEntity)newV);
    PList_delete(delEdge);
    PList_delete(newEdges);
    
    // --- remove regions ---
    temp = NULL;
    while( pRegion region = (pRegion)PList_next(delRegs,&temp) ) {
      M_removeRegion(mesh, region);
    }
    PList_delete(delRegs);
    
    // --- remove and clean faces ---
    temp = NULL;
    while( pFace face = (pFace)PList_next(eFaces,&temp) ) {

      EN_removeData((pEntity)face,"_EdgeSplitMarker_");
      
      for(int iE=0; iE<F_numEdges(face); iE++) {
        pEdge iEdge = F_edge(face,iE);
        if( iEdge != edge ) EN_removeData((pEntity)iEdge,"_EdgeSplitMarker_");
      }

      M_removeFace(mesh, face);
    }
    PList_delete(eFaces);
    
    temp = NULL;
    while( pFace face = (pFace)PList_next(newFaces,&temp) ) {
      EN_removeData((pEntity)face,"_EdgeSplitMarker_");
    }
    PList_delete(newFaces);
    
    // --- remove edge ---
    M_removeEdge(mesh,edge);
    
    return newV ;
  }

  // -------------------------------------------------------------------

}
