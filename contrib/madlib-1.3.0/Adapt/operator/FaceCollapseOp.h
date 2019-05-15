// -*- C++ -*-
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

#ifndef _H_EDGESPLITCOLLAPSEOP
#define _H_EDGESPLITCOLLAPSEOP

#include "MAdOperatorBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class faceCollapseOp: public MAdOperatorBase
  {
  public:

    faceCollapseOp(pMesh, DiscreteSF *);
    faceCollapseOp(const faceCollapseOp &);
    ~faceCollapseOp() {}

    operationType type() const { return MAd_FCOLLAPSE; }

    void collapseOnBoundary(bool, double);
    void reset(pFace, pEdge, bool);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    bool evaluateShapes2D();
    void evaluateLengths() const;

  private:
  
    pFace delFace;
    pEdge delEdge;
    bool clpsOnOppV;
    mutable pVertex oppV;
    double newVPos[3];

    bool constrainBoundary;
    bool checkVol;
    double dVTol; // tolerance on the relative change of volume
  };


  // -------------------------------------------------------------------
  inline faceCollapseOp::faceCollapseOp(const faceCollapseOp & _fc):
    MAdOperatorBase(_fc)
  {
    delFace    = _fc.delFace;
    delEdge    = _fc.delEdge;
    clpsOnOppV = _fc.clpsOnOppV;
    oppV       = _fc.oppV;
    newVPos[0] = _fc.newVPos[0];
    newVPos[1] = _fc.newVPos[1];
    newVPos[2] = _fc.newVPos[2];
    constrainBoundary = _fc.constrainBoundary;
    checkVol = _fc.checkVol;
    dVTol = _fc.dVTol;
  }

  // -------------------------------------------------------------------
  inline faceCollapseOp::faceCollapseOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf)
  {
    delFace = NULL;
    delEdge = NULL;
    clpsOnOppV = true;
    oppV = NULL;
    newVPos[0] = 0.;
    newVPos[1] = 0.;
    newVPos[2] = 0.;
    constrainBoundary = true;
    checkVol = false;
    dVTol = MAdTOL;
  }

  // -------------------------------------------------------------------
  // pFace: face which would disapear
  // pEdge: Edge which would be split
  // clpsOnOppV: Edge collapse type
  //             True if collapse from SplitEdge to existing vertex
  //             False if collapse from existing vertex to SplitEdge
  // -------------------------------------------------------------------
  inline void faceCollapseOp::reset(pFace face, pEdge edge, 
                                    bool typeClps)
  {
    delFace = face;
    delEdge = edge;
    clpsOnOppV = typeClps;

    double Exyz[2][3];
    V_coord(E_vertex(delEdge,0),Exyz[0]);
    V_coord(E_vertex(delEdge,1),Exyz[1]);
  
    for(int i=0; i<3; i++) {
      newVPos[i] = 0.5 * ( Exyz[0][i] + Exyz[1][i] );
    }
  }

  // -------------------------------------------------------------------
  inline void faceCollapseOp::collapseOnBoundary(bool cob, double tolerance)
  {
    constrainBoundary = !cob;
    dVTol = tolerance;
  }
  
  // -------------------------------------------------------------------

}

#endif
