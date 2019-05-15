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

#ifndef _H_FACESWAPOP
#define _H_FACESWAPOP

#include "MAdOperatorBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class faceSwapOp: public MAdOperatorBase
  {
  public:

    faceSwapOp(pMesh, DiscreteSF *);
    faceSwapOp(const faceSwapOp &);
    ~faceSwapOp() {}

    operationType type() const { return MAd_FSWAP; }
  
    void setSwapFace(pFace);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    void evaluateLengths() const;

  private:

    pFace face;           // face to be swapped
    pRegion fRegions[2];  // its regions
    pVertex fOppVerts[2]; // opposite vertices in regions
  };

  // -------------------------------------------------------------------
  inline faceSwapOp::faceSwapOp(const faceSwapOp & _fs):
    MAdOperatorBase(_fs), face(_fs.face)
  {
    fRegions[0]  = _fs.fRegions[0];
    fRegions[1]  = _fs.fRegions[1];
    fOppVerts[0] = _fs.fOppVerts[0];
    fOppVerts[1] = _fs.fOppVerts[1];
  }

  // -------------------------------------------------------------------
  inline faceSwapOp::faceSwapOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf), face(NULL)
  {
    fRegions[0]  = NULL;
    fRegions[1]  = NULL;
    fOppVerts[0] = NULL;
    fOppVerts[1] = NULL;
  }

  // -------------------------------------------------------------------
  inline void faceSwapOp::setSwapFace(pFace _face)
  {
    face = _face;

    for(int iR=0; iR<2; iR++) {
      fRegions[iR] = F_region(face,iR);
      if (fRegions[iR]) fOppVerts[iR] = R_fcOpVt(fRegions[iR],face);
    }

    results->reset();
  }

  // -------------------------------------------------------------------

}

#endif
