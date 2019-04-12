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

#include "VertexMoveOp.h"
#include "OperatorTools.h"
#include "MathUtils.h"

#include <iostream>
using std::cerr;
#include <queue>
using std::queue;
using std::set;
using std::multiset;
#include <cmath>

namespace MAd {

  // -------------------------------------------------------------------
  vDisplacement::vDisplacement(const vDisplacement& _vDisp)
  {
    pv = _vDisp.pv;
    for(int i=0; i<3; i++)  dxyz[i] = _vDisp.dxyz[i];
  }

  // -------------------------------------------------------------------
  vDisplacement::vDisplacement(pVertex v, double disp[3])
  {
    pv = v;
    for(int i=0; i<3; i++)  dxyz[i] = disp[i];
  }

  // -------------------------------------------------------------------
  void vDisplacement::scale(double factor)
  {
    for(int i=0; i<3; i++)  dxyz[i] = factor * dxyz[i];  
  }

  // -------------------------------------------------------------------
  bool vDisplacementLess::operator() (const vDisplacement& vd1,
                                      const vDisplacement& vd2) const
  {
    if (vd1.pv == vd2.pv)
      if (vd1.dxyz[0] == vd2.dxyz[0])
        if (vd1.dxyz[1] == vd2.dxyz[1])
          return (vd1.dxyz[2] < vd2.dxyz[2]);
        else return (vd1.dxyz[1] < vd2.dxyz[1]);
      else return (vd1.dxyz[0] < vd2.dxyz[0]);
    else return (vd1.pv < vd2.pv);

    return false;
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  bool vertexMoveOp::move(std::set<vDisplacement,vDisplacementLess>& _vDisps, 
                          double factor)
  {
    resetDisplacements();
    set<vDisplacement,vDisplacementLess>::const_iterator vDIt  = _vDisps.begin();
    set<vDisplacement,vDisplacementLess>::const_iterator vDEnd = _vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      vDisplacement vDisp = *vDIt;
      vDisp.scale(factor);
      vDisps.insert(vDisp);
    }
    return operate();
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::move(multiset<vDisplacement,vDisplacementLess>& _vDisps, 
                          double factor)
  {
    resetDisplacements();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = _vDisps.begin();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = _vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      vDisplacement vDisp = *vDIt;
      vDisp.scale(factor);
      vDisps.insert(vDisp);
    }
    return operate();
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::move(pVertex v, double disp[3])
  {
    return move( vDisplacement(v,disp) );
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::move(vDisplacement vDisp)
  {
    resetDisplacements();
    addVDisplacement(vDisp);
    return operate();
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::operate()
  {
    double shape;
    if ( !evaluate(&shape) ) return false;
    apply();
    return true;
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::addVDisplacement(pVertex v, double disp[3])
  {
    vDisplacement vDisp(v,disp);
    addVDisplacement(vDisp);
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::addVDisplacement(vDisplacement vDisp)
  {
    vDisps.insert(vDisp);
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::resetDisplacements()
  {
    vDisps.clear();
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::setDisplacement(pVertex v, double disp[3])
  {
    resetDisplacements();
    addVDisplacement(v,disp);
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::setPosition(pVertex v, double pos[3])
  {
    double xyz0[3], dxyz[3];
    V_coord(v,xyz0);
    diffVec(pos,xyz0,dxyz);
    setDisplacement(v,dxyz);
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::setPositionToOptimal(pVertex v)
  {
    double opt[3];
    if (computeOptimalLocation(v, opt)) setPosition(v, opt);
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::computeOptimalLocation(pVertex vt, double opt[3]) const
  {
    if ( V_whatInType(vt) != dim ) return false;

    bool flag = false;

    // --- get original coordinates ---
    double ori[3]; V_coord(vt,ori);

    // --- get connected faces ---
    std::set<pFace> fSet;
    pPList faces = V_faces(vt);
    void * temp = 0;
    while ( pFace pf = (pFace) PList_next(faces,&temp) ) fSet.insert(pf);
    PList_delete(faces);

    // --- get original worst volume ratio ---
    double oriWorstRatio = F_worstVolumeRatio(fSet);
  
    // --- set first trial: the center of the cavity ---
    V_cavityCenter(vt,opt);

    int iter = 0;
    while ( iter < 3 && dotProd(opt,ori) > MAdTOL ) {

      // --- get worst vol ratio if we move to optimal ---
      V_setPosition(vt,opt);
      double optWorstRatio = F_worstVolumeRatio(fSet);

      // --- if it does not improve volume ratio ---
      // --- move it if increasing shape ---
      if( fabs(oriWorstRatio-optWorstRatio) < OPTILOC_MIN_IMPROVE_RATIO ) {

        V_setPosition(vt,ori);
        double oriWorstShp;
        mqm.V_worstShape(vt,&oriWorstShp);

        V_setPosition(vt,opt);
        double optWorstShp;
        mqm.V_worstShape(vt,&optWorstShp);

        if( (optWorstShp-oriWorstShp) > OPTILOC_MIN_IMPROVE_SHAPE ) { 
          flag = true; break;
        }
      }

      // --- if it improves volume ratio ---
      if( (oriWorstRatio-optWorstRatio) > OPTILOC_MIN_IMPROVE_RATIO && 
          optWorstRatio != -1. ) {
        flag = true; break;
      }

      // --- underrelax solution ---
      for(int i=0; i<3; i++)
        opt[i] = 0.5*(opt[i]+ori[i]);

      iter++;
    }

    // --- put back the vertex at its original location ---
    V_setPosition(vt,ori);

    // --- no improvement found ---
    if (!flag) V_coord(vt,opt);

    return flag;
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::checkConstraints() const
  {
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      vDisplacement vd = *vDIt;
      if ( EN_constrained((pEntity)(vd.pv)) )  return false;
    }
    return true;
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::checkGeometry() 
  {
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      vDisplacement vd = *vDIt;
      // check vertex is not on a boundary (if requested)
      if (fixedBndry)    if ( V_whatInType(vd.pv) < dim ) return false;
    }
    return true;
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::evaluateShapes2D()
  {
    pPList faces;
    getCavity(&faces);

    double worstShape = mqm.getElementEvaluator()->bestShapeEver();

    void * temp = NULL;
    while( pFace face = (pFace) PList_next(faces,&temp) ) {

      double fCoords[3][3];
      pMSize fSizes[3] = {NULL,NULL,NULL};

      double fOriNor[3];
      F_normal(face,fOriNor);

      pPList fVerts = F_vertices(face,1);
      void * temp2 = NULL;
      int iNode = 0;
      while ( pVertex pV = (pVertex)PList_next(fVerts,&temp2) ) {

        V_coord(pV,fCoords[iNode]);
      
        multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
        multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
        for (; vDIt != vDEnd; vDIt++) {
          if ( (*vDIt).pv == pV ) {
            for (int iC=0; iC<3; iC++) fCoords[iNode][iC] += (*vDIt).dxyz[iC];
          }
        }
      
        fSizes[iNode] = sizeField->findSize(pV);

        iNode++;
      }
      PList_delete(fVerts);

      // Check if the shape is acceptable and compare it to the worst shape
      double shape = 0.;
      if ( !mqm.getElementEvaluator()->XYZ_F_shape(fCoords, fSizes, fOriNor, &shape) ) {
        PList_delete(faces); 
        return false;
      }
      else {
        if ( shape < worstShape ) worstShape = shape;
      }
    }
    PList_delete(faces);

    results->setWorstShape(worstShape);

    return true;
  }

  // -------------------------------------------------------------------
  bool vertexMoveOp::evaluateShapes()
  {
    if (dim == 2) return evaluateShapes2D();

    pPList regs;
    getCavity(&regs);

    double worstShape = mqm.getElementEvaluator()->bestShapeEver();

    void * temp = NULL;
    while( pRegion region = (pRegion) PList_next(regs,&temp) ) {

      double rCoords[4][3];
      pMSize rSizes[4] = {NULL,NULL,NULL,NULL};

      pPList rVerts = R_vertices(region);
      void * temp2 = NULL;
      int iNode = 0;
      while ( pVertex pV = (pVertex)PList_next(rVerts,&temp2) ) {

        V_coord(pV,rCoords[iNode]);
      
        multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
        multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
        for (; vDIt != vDEnd; vDIt++) {
          if ( (*vDIt).pv == pV ) {
            for (int iC=0; iC<3; iC++) rCoords[iNode][iC] += (*vDIt).dxyz[iC];
          }
        }
      
        rSizes[iNode] = sizeField->findSize(pV);

        iNode++;
      }
      PList_delete(rVerts);

      // Check if the shape is acceptable and compare it to the worst shape
      double shape = 0.;
      if ( !mqm.getElementEvaluator()->XYZ_R_shape(rCoords, rSizes, &shape) ) {
        PList_delete(regs); 
        return false;
      }
      else {
        if ( shape < worstShape ) worstShape = shape;
      }
    }
    PList_delete(regs);

    results->setWorstShape(worstShape);

    return true;
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::evaluateLengths() const
  {
    double minSq = MAdBIG;
    double maxSq = 0.;

    pPList edges;
    getAffectedEdges(&edges);

    void* tmp=0;
    while ( pEdge pe = (pEdge)PList_next(edges,&tmp) ) {
      double lSq = sizeField->SF_E_lengthSq(pe);
      if ( lSq > maxSq ) maxSq = lSq;
      if ( lSq < minSq ) minSq = lSq;
    }

    PList_delete(edges);

    results->setMaxLenSq(maxSq);
    results->setMinLenSq(minSq);
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::getAffectedEdges(pPList * edges) const
  {
    *edges = PList_new();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      pPList vEdges = V_edges(vDIt->pv);
      (*edges) = PList_appPListUnique(*edges,vEdges);
      PList_delete(vEdges);
    }
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::getCavity(pPList * cavity) const
  {
    *cavity = PList_new();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      if ( dim == 3 ) {
        pPList vRegs = V_regions(vDIt->pv);
        (*cavity) = PList_appPListUnique(*cavity,vRegs);
        PList_delete(vRegs);
      }
      else {
        pPList vFaces = V_faces(vDIt->pv);
        (*cavity) = PList_appPListUnique(*cavity,vFaces);
        PList_delete(vFaces);
      }
    }
  }

  // -------------------------------------------------------------------
  void vertexMoveOp::apply()
  {
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDIt  = vDisps.begin();
    multiset<vDisplacement,vDisplacementLess>::const_iterator vDEnd = vDisps.end();
    for (; vDIt != vDEnd; vDIt++) {
      vDisplacement vd = *vDIt;

      double target[3]; V_coord(vd.pv,target);
      for (int i=0; i < 3; i++)  target[i] += vd.dxyz[i];

      if ( !V_setPosition(vd.pv,target) ) {
        cerr << "Error: could not move vertex\n"; throw;
      }
    }
    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}
