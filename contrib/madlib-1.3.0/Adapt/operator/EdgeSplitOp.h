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

#ifndef _H_EDGESPLITOP
#define _H_EDGESPLITOP

#include "MAdOperatorBase.h"
#include "MeshSizeBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class edgeSplitOp: public MAdOperatorBase
  {
  public:

    edgeSplitOp(pMesh, DiscreteSF *);
    edgeSplitOp(const edgeSplitOp &);
    ~edgeSplitOp() { if (xyzSize) delete xyzSize; }

    operationType type() const { return MAd_ESPLIT; }

    // stores the edge to be splitted, computes the location of the split 
    // and returns the (unique) adimensional square length of the resulting edges
    double setSplitEdge(pEdge);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    bool evaluateShapes2D();
    void evaluateLengths() const;

  private:

    pEdge edge;
    double xyz[3], u; // new vertex location: euclidian and edge parameter
    pMSize xyzSize;   // size at the new vertex
  };

  // -------------------------------------------------------------------
  inline edgeSplitOp::edgeSplitOp(const edgeSplitOp & _es):
    MAdOperatorBase(_es)
  {
    edge = _es.edge;
    xyz[0] = _es.xyz[0];
    xyz[1] = _es.xyz[1];
    xyz[2] = _es.xyz[2];
    u = _es.u;
    xyzSize = MS_copy(_es.xyzSize);
  }

  // -------------------------------------------------------------------
  inline edgeSplitOp::edgeSplitOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf), edge(NULL)
  {
    xyz[0] = 0.;
    xyz[1] = 0.;
    xyz[2] = 0.;
    u = -1.;
    xyzSize = NULL;
  }

  // -------------------------------------------------------------------
  inline double edgeSplitOp::setSplitEdge(pEdge _edge)
  {
    if ( xyzSize ) delete xyzSize;
    edge = _edge;
    double lenReduc;
    u = sizeField->SF_E_center(_edge,xyz,&lenReduc,&xyzSize);
    return lenReduc;
  }

  // -------------------------------------------------------------------

}

#endif
