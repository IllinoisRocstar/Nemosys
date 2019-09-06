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

#include "EdgeSplitOp.h"
#include "OperatorTools.h"

namespace MAd {

  // -------------------------------------------------------------------
  bool edgeSplitOp::checkConstraints() const
  {
    if ( EN_constrained((pEntity)edge) )  return false;
    return true;
  }

  // -------------------------------------------------------------------
  bool edgeSplitOp::checkGeometry()
  {
    return true;
  }

  // -------------------------------------------------------------------
  bool edgeSplitOp::evaluateShapes()
  {
    pPList eRegs = E_regions(edge);

    // 2D case
    if ( PList_size(eRegs) == 0 ) {
      bool flag = evaluateShapes2D();
      PList_delete(eRegs);
      return flag;
    }

    // 3D case
    else {

      double worstShape = MAdBIG;

      void * temp = 0;
      while ( pRegion region = (pRegion)PList_next(eRegs,&temp) ) {

        pPList rVerts = R_vertices(region);
      
        // evaluate the two regions resulting from the split
        for( int iR=0; iR<2; iR++ ) {

          pVertex eV = E_vertex(edge,iR);

          pMSize rSizes[4] = { NULL, NULL, NULL, NULL };
          double rCoords[4][3];

          pVertex pV = NULL;
          void * iter = NULL;
          for( int iRV=0; (pV=(pVertex)PList_next(rVerts,&iter) ); iRV++ ) {
            if( pV == eV ) {
              rSizes[iRV] = xyzSize;
              rCoords[iRV][0] = xyz[0]; 
              rCoords[iRV][1] = xyz[1]; 
              rCoords[iRV][2] = xyz[2];
            }
            else {
              rSizes[iRV] = sizeField->findSize(pV);
              V_coord(pV,rCoords[iRV]);
            }
          }
        
          double rShape;
          if( !mqm.getElementEvaluator()->XYZ_R_shape(rCoords,rSizes,&rShape) ) {
            PList_delete (eRegs);
            PList_delete (rVerts);
            return false;
          }

          if( worstShape > rShape ) worstShape = rShape;
        }
        PList_delete(rVerts);
      }

      results->setWorstShape(worstShape);
    }

    PList_delete(eRegs);
  
    return true;
  }

  // -------------------------------------------------------------------
  bool edgeSplitOp::evaluateShapes2D()
  {
    double worstShape = MAdBIG;

    pPList eFaces = E_faces(edge);

    void * temp = NULL;
    while ( pFace face = (pFace)PList_next(eFaces,&temp) ) {

      pPList fVerts = F_vertices(face,1);
    
      // evaluate the two faces resulting from the split
      for( int iF=0; iF<2; iF++ ) {

        pVertex eV = E_vertex(edge,iF);

        pMSize fSizes[3] = { NULL, NULL, NULL };  
        double fCoords[3][3];

        pVertex pV;
        void * iter = NULL;
        for( int iFV=0; (pV=(pVertex)PList_next(fVerts,&iter)); iFV++ ) {
          if( pV == eV ) {
            fSizes[iFV] = xyzSize;
            fCoords[iFV][0] = xyz[0];
            fCoords[iFV][1] = xyz[1];
            fCoords[iFV][2] = xyz[2];
          }
          else {
            fSizes[iFV] = sizeField->findSize(pV);
            V_coord(pV,fCoords[iFV]);
          }
        }
      
        double fShape;
        if( !mqm.getElementEvaluator()->XYZ_F_shape(fCoords,fSizes,0,&fShape) ) {
          PList_delete (eFaces);
          PList_delete (fVerts);
          return false;
        }

        if( worstShape > fShape ) worstShape = fShape;
      }

      PList_delete(fVerts);
    }

    PList_delete(eFaces);

    results->setWorstShape(worstShape);
  
    return true;
  }

  // -------------------------------------------------------------------
  void edgeSplitOp::evaluateLengths() const
  {
    pVertex verts[2];
    double vCoords[2][3];
    pMSize vSizes[2] = {NULL,NULL};
    for( int iV=0; iV<2; iV++ ) {
      verts[iV] = E_vertex(edge,iV);
      V_coord(verts[iV],vCoords[iV]);
      vSizes[iV] = sizeField->findSize(verts[iV]);
    }

    pMSize newSize = sizeField->getSizeOnEntity((pEntity)edge,xyz);

    double lSqMax = sizeField->SF_XYZ_lengthSq(xyz, vCoords[0],
                                               newSize, vSizes[0]);
    double lSqMin = lSqMax;
    double lSq = sizeField->SF_XYZ_lengthSq(xyz, vCoords[1],
                                            newSize, vSizes[1]);  
    if( lSq > lSqMax ) lSqMax = lSq;
    if( lSq < lSqMin ) lSqMin = lSq;
  
    for( int i=0; i<E_numFaces(edge); i++ ) {
      pVertex oppV = F_edOpVt(E_face(edge,i), edge);
      pMSize oppSize = sizeField->findSize(oppV);
      double xyzOpp[3];
      V_coord(oppV,xyzOpp);
      lSq = sizeField->SF_XYZ_lengthSq(xyz, xyzOpp,
                                       newSize, oppSize);
      if( lSq > lSqMax ) lSqMax = lSq;
      if( lSq < lSqMin ) lSqMin = lSq;
    }
  
    if( newSize ) delete newSize;
  
    results->setMaxLenSq(lSqMax);
    results->setMinLenSq(lSqMin);
  }

  // -------------------------------------------------------------------
  void edgeSplitOp::getCavity(pPList * cavity) const
  {
    if ( dim == 3 )  *cavity = E_regions(edge);
    else             *cavity = E_faces(edge);
  }

  // -------------------------------------------------------------------
  void edgeSplitOp::apply()
  {
    E_split(mesh,edge,xyz,u);

    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}
