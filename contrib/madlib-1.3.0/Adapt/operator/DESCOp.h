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

#ifndef _H_DESCOP
#define _H_DESCOP

#include "MAdOperatorBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class DESCOp: public MAdOperatorBase
  {
  public:
  
    DESCOp(pMesh, DiscreteSF *);
    DESCOp(const DESCOp &);
    ~DESCOp() {}

    operationType type() const { return MAd_DESPLTCLPS; }

    void setDESC(pRegion, pEdge, pEdge);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    void evaluateLengths() const;

  private:

    pRegion region;
    pEdge splitE[2];
    double xyz[3];
  };


  // -------------------------------------------------------------------
  inline DESCOp::DESCOp(const DESCOp &_desc):
    MAdOperatorBase(_desc)
  {
    region = _desc.region;
    splitE[0] = _desc.splitE[0];
    splitE[1] = _desc.splitE[1];
  
    xyz[0] = _desc.xyz[0];
    xyz[1] = _desc.xyz[1];
    xyz[2] = _desc.xyz[2];
  }

  // -------------------------------------------------------------------
  inline DESCOp::DESCOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf), region(NULL)
  {
    splitE[0] = NULL;
    splitE[1] = NULL;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
  }

  // -------------------------------------------------------------------

}

#endif
