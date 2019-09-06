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
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_EDGECOLLAPSEOP
#define _H_EDGECOLLAPSEOP

#include "MAdOperatorBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class edgeCollapseOp: public MAdOperatorBase
  {
  public:

    edgeCollapseOp(pMesh, DiscreteSF *);
    edgeCollapseOp(const edgeCollapseOp &);
    ~edgeCollapseOp() {}

    operationType type() const { return MAd_ECOLLAPSE; }

    void collapseOnBoundary(bool, double);
    void setCollapseEdge(pEdge, pVertex, pVertex);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    bool evaluateShapes2D();
    void evaluateLengths() const;
  
  private:

    pEdge edgeDel;
    pVertex vDel; // the vertex to be deleted
    pVertex vTgt; // the target vertex

    bool constrainBoundary;
    double dVTol, dATol; // tolerance on the relative change of volume or area
  };

  // -------------------------------------------------------------------
  inline edgeCollapseOp::edgeCollapseOp(const edgeCollapseOp & _ec):
    MAdOperatorBase(_ec)
  {
    edgeDel = _ec.edgeDel;
    vDel = _ec.vDel;
    vTgt = _ec.vTgt;
    constrainBoundary = _ec.constrainBoundary;
    dVTol = _ec.dVTol;
    dATol = _ec.dATol;
  }

  // -------------------------------------------------------------------
  inline edgeCollapseOp::edgeCollapseOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf)
  {
    edgeDel = NULL;
    vDel = NULL;
    vTgt = NULL;
    constrainBoundary = true;
    dVTol = MAdTOL;
    dATol = MAdTOL;
  }

  // -------------------------------------------------------------------
  inline void edgeCollapseOp::collapseOnBoundary(bool cob, double tolerance)
  {
    constrainBoundary = !cob;
    dVTol = tolerance;
    dATol = tolerance;
  }

  // -------------------------------------------------------------------
  inline void edgeCollapseOp::setCollapseEdge(pEdge _edgeDel,
                                              pVertex _vDel,
                                              pVertex _vTgt)
  {
    edgeDel = _edgeDel;
    vDel = _vDel;
    vTgt = _vTgt;
    results->reset();
  }

  // -------------------------------------------------------------------

}

#endif

