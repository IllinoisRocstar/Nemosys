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

#include "RegionRemoveOp.h"
#include "CallbackManager.h"

namespace MAd {

  // -------------------------------------------------------------------
  bool regionRemoveOp::checkConstraints() const
  {
    return true;
  }

  // -------------------------------------------------------------------
  bool regionRemoveOp::checkGeometry()
  {
    // deny if a vertex is classified on a geometrical edge or vertex
    int nbClassVerts = 0;
    pPList rVerts = R_vertices(region);
    void * temp = NULL;
    while ( pVertex vert = (pVertex) PList_next(rVerts,&temp) ) {
      if ( EN_whatInType((pEntity)vert) < 2 ) nbClassVerts++;
    }
    PList_delete(rVerts);
    if ( nbClassVerts != 0 ) return false;

    // deny if an edge is classified on a geometrical edge
    int nbClassEdges = 0;
    pPList rEdges = R_edges(region);
    void * temp2 = NULL;
    while ( pEdge edge = (pEdge) PList_next(rEdges,&temp2) ) {
      if ( EN_whatInType((pEntity)edge) < 2 ) nbClassEdges++;
    }
    PList_delete(rEdges);
    if ( nbClassEdges != 0 ) return false;

    int nbExternalFaces = 0;
    for (int iF=0; iF < 4; iF++) {

      pFace face = R_face(region,iF);

      // see which faces will be deleted
      if ( F_numRegions(face) == 1 ) {
        nbExternalFaces++;
        delFace[iF] = true;
      }
      else {
        delFace[iF] = false;
      }
    
      // see on which geometric entity the new 
      // boundary faces will be classified
      pGEntity pGE = EN_whatIn((pEntity)face);
      if ( GEN_type(pGE) == 2 ) {
        classifyFaces = true;
        geoFace = pGE;
      }
    }
    // ensure not to create a hole in a volume
    if ( nbExternalFaces == 0 ) return false;

    return false;
  }

  // -------------------------------------------------------------------
  bool regionRemoveOp::evaluateShapes()
  {
    double worstShape = mqm.getElementEvaluator()->bestShapeEver();
    results->setWorstShape( worstShape );
    return true;
  }

  // -------------------------------------------------------------------
  void regionRemoveOp::evaluateLengths() const
  {
    results->setMinLenSq(1.);
    results->setMaxLenSq(1.);
  }

  // -------------------------------------------------------------------
  void regionRemoveOp::apply()
  {
    // --- classify faces ---
    if (classifyFaces) {
      for (int iF=0; iF < 4; iF++) {
        if ( delFace[iF] == false ) {
          pFace face = R_face(region,iF);
          F_setWhatIn(face,geoFace);
        }
      }
    }

    // --- list faces to delete ---
    pPList delF = PList_new();
    for (int iF=0; iF < 4; iF++) {
      if ( delFace[iF] == true ) PList_append(delF,(pEntity)R_face(region,iF));
    }

    // --- list edges to delete ---
    pPList delE = PList_new();
    pPList rEdges = R_edges(region);
    void * temp = NULL;
    while ( pEdge edge = (pEdge) PList_next(rEdges,&temp) ) {
      bool toDel = true;
      pPList eFaces = E_faces(edge);
      void * temp2 = NULL;
      while ( pFace face = (pFace) PList_next(eFaces,&temp2) ) {
        if ( !PList_inList(delF,(pEntity)face) ) { toDel = false; break; }
      }
      PList_delete(eFaces);
      if (toDel) PList_append(delE,(pEntity)edge);
    }
    PList_delete(rEdges);

    // --- list vertices to delete ---
    pPList delV = PList_new();
    pPList rVerts = R_vertices(region);
    void * temp3 = NULL;
    while ( pVertex vert = (pVertex) PList_next(rVerts,&temp3) ) {
      bool toDel = true;
      pPList vEdges = V_edges(vert);
      void * temp4 = NULL;
      while ( pEdge edge = (pEdge) PList_next(vEdges,&temp4) ) {
        if ( !PList_inList(delE,(pEntity)edge) ) { toDel = false; break; }
      }
      if (toDel) PList_append(delV,(pEntity)vert);
    }
    PList_delete(rVerts);
  

    // --- call callback functions ---
    pPList emptyList = PList_new();
    pPList delEn = PList_new();
    void * delTmp = NULL;
    while ( pFace face = (pFace) PList_next(delF,&delTmp) ) PList_append(delEn,(pEntity)face);
    delTmp = NULL;
    while ( pEdge edge = (pEdge) PList_next(delE,&delTmp) ) PList_append(delEn,(pEntity)edge);
    delTmp = NULL;
    while ( pVertex vert = (pVertex) PList_next(delV,&delTmp) ) PList_append(delEn,(pEntity)vert);
    CallBackManagerSgl::instance().callCallBacks(delEn,
                                                 emptyList,
                                                 MAd_RREMOVE,
                                                 (pEntity)region);
    PList_delete(delEn);
  
    // --- delete entities ---

    M_removeRegion(mesh,region);

    void * delTmp2 = NULL;
    while ( pFace face = (pFace) PList_next(delF,&delTmp2) ) {
      M_removeFace(mesh,face);
    }

    delTmp2 = NULL;
    while ( pEdge edge = (pEdge) PList_next(delE,&delTmp2) ) {
      M_removeEdge(mesh,edge);
    }

    delTmp2 = NULL;
    while ( pVertex vertex = (pVertex) PList_next(delV,&delTmp2) ) {
      M_removeVertex(mesh,vertex);
    }

    PList_delete(delF);
    PList_delete(delE);
    PList_delete(delV);

    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}
