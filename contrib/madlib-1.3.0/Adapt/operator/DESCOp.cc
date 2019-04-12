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

#include "DESCOp.h"
#include "OperatorTools.h"
#include "MeshSizeBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  void DESCOp::setDESC(pRegion r, pEdge e0, pEdge e1)
  {
    region = r;
    splitE[0] = e0;
    splitE[1] = e1;
    E_center(splitE[0],xyz);  // Collapse on mid-edge e0 only 
  }

  // -------------------------------------------------------------------
  bool DESCOp::checkConstraints() const
  {
    if ( EN_constrained((pEntity)splitE[0]) || 
         EN_constrained((pEntity)splitE[1]) ) {
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------
  bool DESCOp::checkGeometry()
  {
    if( E_whatInType(splitE[1]) != 3 ) return false;
  
    for( int i=0; i<2; i++ ) {

      int nbV = 0;
      pVertex verts[2];
      for(int iF=0; iF < R_numFaces(region); iF++) {
        pFace face = R_face(region,iF);
        if( F_inClosure(face,(pEntity)splitE[i]) ) {
          pRegion oppR = F_region(face,0);
          if( oppR == region ) oppR = F_region(face,1);
          if ( oppR ) {
            verts[nbV] = R_fcOpVt(oppR,face);
            nbV++;
          }
        }
      }

      pEdge edge = splitE[(i+1)%2];
      for( int j=0; j<E_numFaces(edge); j++ ) {
        pVertex oppV = F_edOpVt(E_face(edge,j),edge);
        for( int k=0; k<nbV; k++ ) {
          if( oppV == verts[k] ) return false;
        }
      }
    }

    return true;
  }

  // -------------------------------------------------------------------
  bool DESCOp::evaluateShapes()
  {
    double worstShape = mqm.getElementEvaluator()->bestShapeEver();

    pMSize xyzSize = sizeField->getSizeOnEntity((pEntity)region,xyz);

    for(int iE=0; iE<2; iE++) {
    
      pPList eRegs = E_regions(splitE[iE]);

      void * temp = NULL;
      while( pRegion pR = (pRegion)PList_next(eRegs,&temp) ) {
    	
        if( pR == region ) continue;
	
        pPList rVerts = R_vertices(pR);

        // evaluate the two regions resulting from the split
        for( int iER=0; iER<2; iER++ ) {

          pVertex eV = E_vertex(splitE[iE],iER);

          pMSize rSizes[4];  
          for (int idel = 0; idel<4; idel++) rSizes[idel] = NULL;
          double rCoords[4][3];

          pVertex pV = NULL;
          void * iter = NULL;
          for( int iERV=0; ( pV=(pVertex)PList_next(rVerts,&iter) ); iERV++ ) {
            if( pV == eV ) {
              rSizes[iERV] = xyzSize;
              rCoords[iERV][0] = xyz[0]; 
              rCoords[iERV][1] = xyz[1]; 
              rCoords[iERV][2] = xyz[2];
            }
            else {
              rSizes[iERV] = sizeField->findSize(pV);
              V_coord(pV,rCoords[iERV]);
            }
          }
        
          double rShape;
          if( !mqm.getElementEvaluator()->XYZ_R_shape(rCoords,rSizes,&rShape) ) {
            PList_delete (rVerts);
            if( xyzSize ) delete xyzSize;
            PList_delete(eRegs);
            return false;
          }
	  
          if( worstShape > rShape ) worstShape = rShape;
        }

        PList_delete(rVerts);
      }

      PList_delete(eRegs);
    }
  
    if( xyzSize ) delete xyzSize;

    results->setWorstShape(worstShape); 

    return true;
  }

  // -------------------------------------------------------------------
  void DESCOp::evaluateLengths() const
  {
    // get size at new point (interpolation)
    pMSize xyzSize = sizeField->getSizeOnEntity((pEntity)splitE[0],xyz);

    double lSqMin = MAdBIG, lSqMax = 0.;

    for( int iE=0; iE<2; iE++ ) {
      for( int i=0; i<E_numFaces(splitE[iE]); i++ ) {

        pVertex oppV = F_edOpVt( E_face(splitE[iE],i), splitE[iE] );
        const pMSize oppSize = sizeField->findSize(oppV);
      
        double xyzOpp[3];
        V_coord(oppV,xyzOpp);
        double lSq = sizeField->SF_XYZ_lengthSq(xyz, xyzOpp,
                                                xyzSize, oppSize);
        if( lSq > lSqMax ) lSqMax = lSq;
        if( lSq < lSqMin ) lSqMin = lSq;
      }
    }

    results->setMaxLenSq(lSqMax);
    results->setMinLenSq(lSqMin);
  
    if( xyzSize ) delete xyzSize;
  }


  // -------------------------------------------------------------------
  void DESCOp::getCavity(pPList * cavity) const
  {
    *cavity = PList_new();
    PList_append(*cavity,(pEntity)region);

    for(int i=0; i<2; i++) {
      pPList eRegs = E_regions(splitE[i]);
      void * iter = NULL;
      for(pRegion pR; ( pR=(pRegion)PList_next(eRegs,&iter) ); ) {
        if( pR != region ) {
          PList_append(*cavity,(pEntity)pR);
        }
      }
      PList_delete(eRegs);
    }
  }

  // -------------------------------------------------------------------
  void DESCOp::apply()
  {
    pVertex newV[2];  
    for(int iE=0; iE<2; iE++) {
      newV[iE] = E_split(mesh, splitE[iE], xyz);
    }

    pVertex vTgt = newV[0];
  
    pEdge edgeDel = E_exist(vTgt,newV[1]);
    assert ( edgeDel );
  
    E_collapse(mesh, edgeDel, newV[1], vTgt);
 
    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}
