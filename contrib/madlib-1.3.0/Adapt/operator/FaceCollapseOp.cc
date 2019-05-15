// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Olivier Pierard, Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#include "FaceCollapseOp.h"
#include "OperatorTools.h"
#include "MeshSizeBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  bool faceCollapseOp::checkConstraints() const
  {
    // find the opposite vertex
    oppV = F_edOpVt(delFace,delEdge);

    if( EN_constrained((pEntity)delEdge)) return false;
    if( !clpsOnOppV && EN_constrained((pEntity)oppV) ) return false;

    return true;
  }

  // -------------------------------------------------------------------
  bool faceCollapseOp::checkGeometry()
  {
    pGEntity fGE = F_whatIn(delFace);
    if( clpsOnOppV )
      {
        // delEdge and delFace must have the same classification
        if( E_whatIn(delEdge) != fGE ) return false;
      }
    else
      {
        // oppV and delFace must have the same classification
        if( V_whatIn(oppV) != fGE ) return false;
      }
    
    // It leaves us with two possibilities:
    // (a) In case ( clpsOnOppV == true ) :
    //     (1) delFaces and delEdges are classified on a dim 2
    //     (2) delFaces and delEdges are classified on a dim 3
    // (b) In case ( clpsOnOppV == false ) :
    //     (1) delFaces and oppV are classified on a dim 2
    //     (2) delFaces and oppV are classified on a dim 3
    // For both (a) and (b):
    // In case (2) no dimension reduction can occur -> unconditionally valid
    // In case (1): * if 2D mesh, no dim reduc -> unconditionally valid
    //              * if 3D mesh, two things to check:
    //                   A.: check boundaries are not constrained
    //                       (otherwise ok but checkVol=true)
    //                   B.: check dimension reduction

    if ( dim == 3 && GEN_type(fGE) == 2 )
      {
        // check A
        if( constrainBoundary ) return false;
        else checkVol = true;
            
        // check B
        // - Dimension reduction on node (on oppV, case clpsOnOppV=true): as delEdge 
        //   is classified on dim 2, only one surface (fGE) passes 
        //   through it. As it is the same as the geom ent of delFace, 
        //   oppV already touches it -> no node dim reduction.
        // - Dimension reduction on lines and surfaces:
        //   occurs if these 2 conditions hold together:
        //   C1: for one of the regions bounding delFace, the face 
        //       containing delEdge and which is not delFace is 
        //       classified on a surface,
        //   C2: for one of the regions bounding delFace, the edge 
        //       of the region opposite to delEdge is classified 
        //       on a dimension < 3 (line contact).
        //   C2b (surface contact): for one of the regions bounding 
        //       delFace, one of the 2 other faces is classified on 
        //       a surface. This is found by Condition 2.
        // - Dimension reduction occur if these 2 conditions hold together:
        //   C3: the face has only 1 region,
        //   C4: there is a face using delEdge which opposite vertex
        //       is linked to oppV by an edge oppE and oppE is not 
        //       used by the only region.
        bool C1 = false, C2 = false;
        for (int iR=0; iR<2; iR++)
          {
            pRegion reg = F_region(delFace,iR);
            if ( !reg ) continue;
            
            // check C1
            for (int iF=0; iF<4; iF++) {
              pFace pf = R_face(reg,iF);
              if ( pf == delFace ) continue;
              if ( F_inClosure(pf,(pEntity)delEdge) && F_whatInType(pf) == 2 ) C1 = true;
            }

            // check C2
            if ( E_whatInType( R_gtOppEdg(reg, delEdge) ) < 3 ) C2 = true;
          }
        if ( C1 && C2 ) return false;

        // In case C3 holds
        if ( F_numRegions(delFace) == 1 )
          {
            // check C4
            pRegion pr = F_region(delFace,0);
            if ( !pr ) pr = F_region(delFace,1);
            pPList delEdgeFaces = E_faces(delEdge);
            void * temp = NULL;
            pFace pf;
            while ( ( pf = (pFace) PList_next(delEdgeFaces,&temp) ) )
              {
                if ( pf == delFace ) continue;
                if ( R_inClosure(pr,(pEntity)pf) ) continue;
                if ( E_exist(oppV, F_edOpVt(pf,delEdge) ) ) {
                  PList_delete(delEdgeFaces); return false;
                }
              }
            PList_delete(delEdgeFaces);
          }
      }

    return true;
  }

  // -------------------------------------------------------------------
  bool faceCollapseOp::evaluateShapes()
  {
    if( dim != 3 ) return evaluateShapes2D();
  
    double worstShape = mqm.getElementEvaluator()->bestShapeEver();

    pRegion delReg[2];
    delReg[0] = F_region(delFace,0);
    delReg[1] = F_region(delFace,1);

    double rCoords[4][3];   
    pMSize rSizes[4] = {NULL, NULL, NULL, NULL};

    double volIn = 0.;
    double volOut = 0.;
  
    // 1°) General case valid for any clpsOnOppV
  
    pPList eRegs = E_regions(delEdge);
    void * temp = NULL;
    while( pRegion region = (pRegion)PList_next( eRegs, &temp ) ) {   
    
      if (checkVol) { volIn += R_volume (region); }
    
      if( region==delReg[0] || region==delReg[1] ) continue;
    
      // check every tet resulting from the split
      for( int iR=0; iR<2; iR++ ) {
      
        pVertex eV = E_vertex(delEdge,iR);
      
        // Loop over all the vertices of region
        pVertex rV;
        for( int iRV=0; iRV < R_numVertices(region); iRV++ ) {        
        
          rV = R_vertex( region, iRV );
      
          if ( rSizes[iRV] ) delete rSizes[iRV];
        
          if( rV == eV ) { // Set oppV size and pos. instead of end of delEdge
            if( clpsOnOppV ) {  // Collapse on oppV
              double xyz[3]; 
              V_coord(oppV,xyz);
              for( int k=0; k<3; k++ ) rCoords[iRV][k] = xyz[k];
              rSizes[iRV] = sizeField->getSize(oppV);  
            }
            else {  // Collapse on newVt
              for( int k=0; k<3; k++ ) rCoords[iRV][k] = newVPos[k];
              rSizes[iRV] = sizeField->getSizeOnEntity((pEntity)delEdge,rCoords[iRV]);
            }
          }
          else {
            rSizes[iRV] = sizeField->getSize(rV);
            V_coord(rV,rCoords[iRV]);
          }
        }
      
        // Some cleaning and exit if new element is not acceptable
        double shape;
        if( ! mqm.getElementEvaluator()->XYZ_R_shape(rCoords,rSizes,&shape) ) {
          for (int i=0; i<4; i++) {
            if ( rSizes[i] ) delete rSizes[i];
          }
          PList_delete(eRegs);
          return false;
        }

        if( worstShape > shape ) worstShape = shape;
      
        if (checkVol) { volOut += R_XYZ_volume(rCoords); }
      }  // End of loop on sub-tetra
    }  // End loop on regions around delEdge
    PList_delete(eRegs);
  
    // 2°) Particular case if collapse edge on new vertex; need other checks 
  
    if( !clpsOnOppV ) {     
    
      // loop over the regions of oppV excepted delReg[0-1]
      pPList vRegs = V_regions(oppV);
      void * temp = NULL;
      while( pRegion region = (pRegion)PList_next(vRegs,&temp) ) {
    
        if( region==delReg[0] || region==delReg[1] ) continue;
      
        if (checkVol) { volIn += R_volume(region); }
      
        pVertex rV;
        for(int iRV=0; iRV < R_numVertices(region); iRV++) {
      
          rV = R_vertex(region, iRV);
      
          if ( rSizes[iRV] ) delete rSizes[iRV];
        
          if( rV == oppV ) { // Set size and pos. of new vertex on delEdge
            for(int iC=0; iC<3; iC++) rCoords[iRV][iC] = newVPos[iC];
            rSizes[iRV] = sizeField->getSizeOnEntity((pEntity)delEdge,rCoords[iRV]);
          }
          else {
            rSizes[iRV] = sizeField->getSize(rV);
            V_coord(rV,rCoords[iRV]);
          }
        }
      
        if (checkVol) { volOut += R_XYZ_volume(rCoords); }
      
        // Some cleaning and exit if new element is not acceptable
        double shape;
        if( ! mqm.getElementEvaluator()->XYZ_R_shape(rCoords,rSizes,&shape) ) {
          for (int i=0; i<4; i++) {
            if ( rSizes[i] ) delete rSizes[i];
          }
          PList_delete(vRegs);
          return false;
        }
      
        if( worstShape > shape ) worstShape = shape;
      
      }  // End loop on vRegs
      PList_delete(vRegs);
    }
  
    for (int i=0; i<4; i++) {
      if ( rSizes[i] ) delete rSizes[i];
    }
  
    results->setWorstShape(worstShape); 

    if (checkVol) {
      double dV = fabs( (volIn-volOut) / volIn );
      if( dV > dVTol ) return false;
    }
  
    return true;
  }

  // -------------------------------------------------------------------
  bool faceCollapseOp::evaluateShapes2D()
  {
    double worstShape = mqm.getElementEvaluator()->bestShapeEver();

    // get the normal to the del face
    double fNormal[3];
    F_normal(delFace,fNormal);

    double fCoords[3][3];   
    pMSize fSizes[3] = {NULL, NULL, NULL};
  
    // 1°) General case valid for any clpsOnOppV
  
    for( int eF=0; eF < E_numFaces(delEdge); eF++ ) {
    
      pFace face = E_face(delEdge, eF);
    
      if( face == delFace ) continue;

      for( int iSubF=0; iSubF<2; iSubF++ ) {
      
        pVertex eV = E_vertex(delEdge,iSubF);
      
        // Loop over all the vertices of face
        pVertex fV;
        for( int iFV=0; iFV<F_numVertices(face); iFV++ ) {        
        
          fV = F_vertex(face, iFV);
      
          if ( fSizes[iFV] ) delete fSizes[iFV];
        
          if( fV == eV ) { // Set oppV size and pos. instead of end of delEdge
            if( clpsOnOppV ) {
              double xyz[3]; 
              V_coord(oppV,xyz);
              for(int k=0; k<3; k++ ) fCoords[iFV][k] = xyz[k];
              fSizes[iFV] = sizeField->getSize(oppV);  
            }
            else {  // Collapse on newVt
              for(int k=0; k<3; k++) fCoords[iFV][k] = newVPos[k];
              fSizes[iFV] = sizeField->getSizeOnEntity((pEntity)delEdge,fCoords[iFV]);
            }
          }
          else {
            V_coord(fV,fCoords[iFV]);
            fSizes[iFV] = sizeField->getSize(fV);
          }
        }
      
        // Some cleaning and exit if new element is not acceptable
        double shape;
        if( ! mqm.getElementEvaluator()->XYZ_F_shape(fCoords,fSizes,fNormal,&shape) ) {
          for (int i = 0; i<3; i++) {
            if ( fSizes[i] ) delete fSizes[i];
          }
          return false;
        }
      
        if( worstShape > shape ) worstShape = shape;
      }  // End of loop on sub-tetra
    }  // End loop on regions around delEdge
  
  
    // 2°) Particular case if collapse edge on new vertex; need other checks 
  
    if( !clpsOnOppV ) {
    
      pPList vFaces = V_faces(oppV);
      void * temp = NULL;
      while( pFace face = (pFace)PList_next( vFaces, &temp ) ) { 
    
        if( face == delFace ) continue;
      
        pVertex fV;       
        for( int k=0; k < F_numVertices(face); k++ ) {
      
          fV = F_vertex(face,k);
      
          if ( fSizes[k] ) delete fSizes[k];
        
          if( fV == oppV ) { // Set size and pos. of new vertex on delEdge
            for( int iC=0; iC<3; iC++ ) fCoords[k][iC] = newVPos[iC];
            fSizes[k] = sizeField->getSizeOnEntity((pEntity)delEdge,fCoords[k]);
          }
          else {
            V_coord(fV,fCoords[k]);
            fSizes[k] = sizeField->getSize(fV);
          }
        }
      
        // Some cleaning and exit if new element is not acceptable
        double shape;
        if( ! mqm.getElementEvaluator()->XYZ_F_shape(fCoords,fSizes,fNormal,&shape) ) {
          for (int i = 0; i<3; i++) {
            if ( fSizes[i] ) delete fSizes[i];
          }
          return false;
        }
      
        if( worstShape > shape ) worstShape = shape;
      }  // End loop on vFaces              
      PList_delete(vFaces);
    }
  
    for (int i = 0; i<3; i++) {
      if ( fSizes[i] ) delete fSizes[i];
    }
  
    results->setWorstShape(worstShape);
    
    return true;
  }

  // -------------------------------------------------------------------
  void faceCollapseOp::evaluateLengths() const
  {
    double minSq = MAdBIG;
    double maxSq = 0.;
  
    // Interpolate the desired size at splitting location (for !clpsOnOppV)
    pMSize sizeNewV = NULL;
    if (!clpsOnOppV) {
      sizeNewV = sizeField->getSizeOnEntity((pEntity)delEdge,newVPos);
    }
  
    // 1°) General case, all faces touching delEdge are modified for both types
  
    pFace face;
    for( int i=0; i<E_numFaces(delEdge); i++ ) {
    
      face = E_face(delEdge, i);
    
      if( face == delFace ) continue;
    
      pVertex vertex = F_edOpVt(face,delEdge);
    
      double lenSq;
      if( clpsOnOppV ) {
        lenSq = sizeField->SF_VV_lengthSq(oppV,vertex);
      }
      else {
        double vxyz[3];
        V_coord(vertex,vxyz);
        pMSize vSize = sizeField->findSize(vertex);
        lenSq = sizeField->SF_XYZ_lengthSq( newVPos, vxyz, 
                                            sizeNewV, vSize );
      }

      if( lenSq > maxSq ) maxSq = lenSq;
      if( lenSq < minSq ) minSq = lenSq;
    }
    
    // 2°) For type 0 only, all edges touching oppV are affected as well
    
    if( !clpsOnOppV ) {
    
      for( int j=0; j<V_numEdges( oppV ); j++ )
        {
          pVertex vertex = E_otherVertex( V_edge(oppV, j), oppV );
          double vxyz[3];
          V_coord(vertex,vxyz);
          pMSize vSize = sizeField->findSize(vertex);
          double lenSq = sizeField->SF_XYZ_lengthSq( newVPos, vxyz, 
                                                     sizeNewV, vSize );
      
          if( lenSq > maxSq ) maxSq = lenSq;
          if( lenSq < minSq ) minSq = lenSq;
        }          
    }
  
    results->setMaxLenSq(maxSq);
    results->setMinLenSq(minSq);

    if (sizeNewV) delete sizeNewV;
  }

  // -------------------------------------------------------------------
  void faceCollapseOp::getCavity(pPList * cavity) const
  {
    if( clpsOnOppV )  // Delete new vertex
      {
        if ( dim == 3 ) *cavity = E_regions(delEdge);
        else            *cavity = E_faces(delEdge);
      }
    else  // Delete existing vertex
      {
        if ( dim == 3 ) {
          *cavity = V_regions(oppV);
          pPList eRegs = E_regions(delEdge);
          void * temp = 0;
          while ( pRegion region = (pRegion)PList_next(eRegs,&temp) ) {
            if ( !PList_inList(*cavity,(pEntity)region) ) PList_append(*cavity,(pEntity)region);
          }
          PList_delete(eRegs);
        }
        else {
          *cavity = V_faces(oppV);
          pPList eFaces = E_faces(delEdge);
          void * temp = 0;
          while ( pFace face = (pFace)PList_next(eFaces,&temp) ) {
            if ( !PList_inList(*cavity,(pEntity)face) ) PList_append(*cavity,(pEntity)face);
          }
          PList_delete(eFaces);

        }
      }
  }

  // -------------------------------------------------------------------
  void faceCollapseOp::apply()
  {
    pVertex newV = E_split(mesh, delEdge, newVPos);
  


    // determine the edge collapse
    pEdge edgeDel = E_exist(oppV,newV);
    pVertex vDel, vTgt;
    if( clpsOnOppV ) { vDel = newV; vTgt = oppV; } 
    else             { vDel = oppV; vTgt = newV; }

    E_collapse(mesh,edgeDel,vDel,vTgt);

    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}
