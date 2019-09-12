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

#ifndef _H_VERTEXMOVEOP
#define _H_VERTEXMOVEOP

#include "MAdOperatorBase.h"
#include "MobileObject.h"
#include "DistanceFunction.h"

#include <set>

#define OPTILOC_MIN_IMPROVE_RATIO  1.e-3
#define OPTILOC_MIN_IMPROVE_SHAPE  1.e-3

namespace MAd {

  // -------------------------------------------------------------------
  class vertexMoveOp : public MAdOperatorBase
  {
  public:

    vertexMoveOp(pMesh, DiscreteSF *, bool _fixBndry=true);
    vertexMoveOp(const vertexMoveOp &);
    ~vertexMoveOp() {}

    operationType type() const { return MAd_VERTEXMOVE; }

    void setFixedBndry(bool);

    void resetDisplacements();
    void addVDisplacement(pVertex, double[3]);
    void addVDisplacement(vDisplacement);
    void setDisplacement(pVertex, double[3]);
    void setPosition(pVertex, double[3]);
    void setPositionToOptimal(pVertex);

    bool computeOptimalLocation(pVertex, double[3]) const;

    void getAffectedEdges(pPList *) const;
    void getCavity(pPList *) const;

    void apply();
    bool move (std::set<vDisplacement,vDisplacementLess>&, double factor=1.0);
    bool move (std::multiset<vDisplacement,vDisplacementLess>&, double factor=1.0);
    bool move (pVertex, double[3]);
    bool move (vDisplacement);
    bool operate ();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    bool evaluateShapes2D();
    void evaluateLengths() const;

  private:

    std::multiset<vDisplacement,vDisplacementLess> vDisps;
    bool fixedBndry; // indicates if the boundary vertices can move
  };

  // -------------------------------------------------------------------
  inline vertexMoveOp::vertexMoveOp(const vertexMoveOp & _vm):
    MAdOperatorBase(_vm) 
  {
    vDisps = _vm.vDisps;
    fixedBndry = _vm.fixedBndry;
  }

  // -------------------------------------------------------------------
  inline vertexMoveOp::vertexMoveOp(pMesh _m, DiscreteSF * _sf, 
                                    bool _fixBndry): 
    MAdOperatorBase(_m,_sf), fixedBndry(_fixBndry)
  {}

  // -------------------------------------------------------------------
  inline void vertexMoveOp::setFixedBndry(bool fixed)
  {
    fixedBndry = fixed;
  }

  // -------------------------------------------------------------------

}

#endif
