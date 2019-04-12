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
// Authors: Arnaud Francois, Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_EDGESWAPOP
#define _H_EDGESWAPOP

#include "MAdOperatorBase.h"
#include "EdgeSwapConfig.h"

namespace MAd {

  // -------------------------------------------------------------------
  class edgeSwapOp : public MAdOperatorBase
  {
  public:

    edgeSwapOp();
    edgeSwapOp(pMesh, DiscreteSF *);
    edgeSwapOp(const edgeSwapOp &);
    ~edgeSwapOp() {}

    operationType type() const { return MAd_ESWAP; }
    int getMaxNumRgns() { return 7; }

    void setSwapEdge(pEdge);
    void swapOnBoundary(bool,double);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    void evaluateLengths() const;

    bool checkGeometry2D();
    bool checkGeometry3D();
    bool evaluateShapes2D();
    bool evaluateShapes3D();

    int apply3D();
    int apply2D();

    pEdge edge;          // the edge to be swapped
    int conf;            // best swap configuration - set by geomCheck()
    pVertex vt;          // top and bottom vertices of the edge (define the polygon orientation
    pVertex vb;          // w/r to vertices)
    std::vector<pVertex> vertices; // oriented list of vertices defining the polygon

    EdgeSwapConfiguration swap_config;

    bool constrainBoundary;
    bool checkVol; // tells if it's interesting to compute the variation of 
                   // volume in the cavity (deduced from classif of edge 
                   // and constrainBoundary)
    double dVTol;

    // reporting
    void reportSwap( );
    bool   report;
    int    reportId;
    std::string reportPrefix;

  };

  // -------------------------------------------------------------------
  inline edgeSwapOp::edgeSwapOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf), edge(0), conf(-1), vt(0), vb(0), 
    constrainBoundary(true),  checkVol(false), dVTol(MAdTOL),
    report(false), reportId(0), reportPrefix("")
  {
    vertices.reserve(7);
  }

  // -------------------------------------------------------------------
  inline edgeSwapOp::edgeSwapOp(const edgeSwapOp & _es):
    MAdOperatorBase(_es), swap_config( _es.swap_config )
  {
    edge        = _es.edge;
    conf        = _es.conf;
    vt          = _es.vt;
    vb          = _es.vb;
    vertices    = _es.vertices;

    constrainBoundary = _es.constrainBoundary;
    checkVol          = _es.checkVol;
    dVTol             = _es.dVTol;

    report       = _es.report;
    reportId     = _es.reportId;
    reportPrefix = _es.reportPrefix;
  }

  // -------------------------------------------------------------------
  inline void edgeSwapOp::setSwapEdge(pEdge e) 
  {
    edge = e;
    results->reset();
    conf = -1;
    checkVol = false;
  }
  
  // -------------------------------------------------------------------
  inline void edgeSwapOp::swapOnBoundary(bool eswapob, double dV)
  {
    constrainBoundary = !eswapob;
    dVTol = dV;
  }
}
// -------------------------------------------------------------------
#endif
