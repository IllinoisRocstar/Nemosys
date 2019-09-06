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

#include "FaceSwapOp.h"
#include "CallbackManager.h"

namespace MAd {

  // -------------------------------------------------------------------
  bool faceSwapOp::checkConstraints() const
  {
    if ( EN_constrained((pEntity)face) ) return false;
    return true;
  }

  // -------------------------------------------------------------------
  bool faceSwapOp::checkGeometry()
  {
    if ( F_whatInType(face) != 3 )       return false;
    if ( !fRegions[0] || !fRegions[1] )  return false;
    if ( R_whatIn(fRegions[0]) != R_whatIn(fRegions[1]) ) return false;
    return true;
  }

  // -------------------------------------------------------------------
  bool faceSwapOp::evaluateShapes()
  {
    // this case would lead to an edge swap
    if ( E_exist(fOppVerts[0],fOppVerts[1]) ) return false;

    double worst = mqm.getElementEvaluator()->bestShapeEver();

    pPList r0Verts = R_vertices(fRegions[0]);

    // collect coordinates of the five vertices
    double xyzR0[4][3];
    void * temp1 = 0; 
    int iV = 0;
    while( pVertex pV = (pVertex)PList_next(r0Verts,&temp1) ) {
      V_coord(pV,xyzR0[iV]);
      iV++;
    }
    double xyzV1[3];
    V_coord(fOppVerts[1],xyzV1);

    // calculate worst shape of the three new regions
    pPList fVerts = F_vertices(face,1);
    void * temp2 = 0;
    while( pVertex fVertex = (pVertex)PList_next(fVerts,&temp2) ) {

      double xyzNewR[4][3];
      pMSize newRSizes[4] = { NULL, NULL, NULL, NULL };
    
      void * temp3 = 0;
      int iR0V = 0;
      while( pVertex pVR0 = (pVertex)PList_next(r0Verts,&temp3) ) {
        if(pVR0 == fVertex) { 
          newRSizes[iR0V] = sizeField->findSize(fOppVerts[1]);
          for (int iC=0; iC<3; iC++) xyzNewR[iR0V][iC] = xyzV1[iC];
          iR0V++;
        }
        else {
          newRSizes[iR0V] = sizeField->findSize(pVR0);
          for (int iC=0; iC<3; iC++) xyzNewR[iR0V][iC] = xyzR0[iR0V][iC];
          iR0V++;
        }
      }
    
      double shape;
      if( !mqm.getElementEvaluator()->XYZ_R_shape(xyzNewR,newRSizes,&shape) ) {
        PList_delete(fVerts);
        PList_delete(r0Verts);
        return false;
      }
      if(shape < worst) worst = shape;
    }
    PList_delete(fVerts);
    PList_delete(r0Verts);
  
    results->setWorstShape(worst);
  
    return true;
  }

  // -------------------------------------------------------------------
  void faceSwapOp::evaluateLengths() const
  {
    double xyz[2][3];
    pMSize sizes[2] = { NULL, NULL };
    for(int iV=0; iV<2; iV++) {
      V_coord(fOppVerts[iV],xyz[iV]);
      sizes[iV] = sizeField->findSize(fOppVerts[iV]);
    }

    double min = sizeField->SF_XYZ_lengthSq(xyz[0], xyz[1],
                                            sizes[0], sizes[1]);
    double max = min;

    results->setMinLenSq(min);
    results->setMaxLenSq(max);
  }

  // -------------------------------------------------------------------
  void faceSwapOp::getCavity(pPList * cavity) const
  {
    *cavity = PList_new();
    PList_append(*cavity, (pEntity)(fRegions[0]));
    PList_append(*cavity, (pEntity)(fRegions[1]));
  }

  // -------------------------------------------------------------------
  void faceSwapOp::apply()
  {
    pPList newRegions = PList_new();

    // get the geometric entity on which the cavity is classified
    pGRegion geoEntity = R_whatIn(fRegions[0]);

    // create new edge
    pEdge newEdge = M_createE(mesh,fOppVerts[0],fOppVerts[1],(pGEntity)geoEntity);

    // create new faces
    pFace newFaces[3];
    pPList fVerts = F_vertices(face,1);
    void * temp = 0;
    int iF = 0;
    while ( pVertex pV = (pVertex)PList_next(fVerts,&temp) ) {
    
      pEdge fEdges[3];
      fEdges[0] = E_exist(fOppVerts[0],pV);
      fEdges[1] = E_exist(pV,fOppVerts[1]);
      fEdges[2] = newEdge;

      newFaces[iF] = M_createF(mesh,3,fEdges,(pGEntity)geoEntity);

      iF++;
    }

    pFace boundaryFaces[2][3];
    for(int iR=0; iR<2; iR++) {
      for(int iV=0; iV<3; iV++) {
        boundaryFaces[iR][iV] = R_vtOpFc(fRegions[iR],(pVertex)PList_item(fVerts,iV));
      }
    }
    PList_delete(fVerts);

    // create new regions
    for (int iR=0; iR<3; iR++) {

      pFace rFaces[4];
      rFaces[0] = boundaryFaces[0][iR];
      rFaces[1] = boundaryFaces[1][iR];
      rFaces[2] = newFaces[(iR+1)%3];
      rFaces[3] = newFaces[(iR+2)%3];

      pRegion newR = M_createR(mesh,4,rFaces,(pGEntity)geoEntity);
      PList_append(newRegions, (pEntity)newR);
    }

    // list deleted regions
    pPList oldRegions = PList_new();
    PList_append(oldRegions,(pEntity)(fRegions[0]));
    PList_append(oldRegions,(pEntity)(fRegions[1]));

    // call callback functions
    CallBackManagerSgl::instance().callCallBacks(oldRegions,
                                                 newRegions,
                                                 MAd_FSWAP,
                                                 (pEntity)face);
    PList_delete(oldRegions);
    PList_delete(newRegions);

    // delete old cavity
    M_removeRegion(mesh,fRegions[0]);
    M_removeRegion(mesh,fRegions[1]);
    M_removeFace(mesh,face);

    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}


