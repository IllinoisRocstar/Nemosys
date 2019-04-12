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

#ifndef _H_REGIONREMOVEOP
#define _H_REGIONREMOVEOP

#include "MeshDataBaseInterface.h"
#include "MAdOperatorBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class regionRemoveOp: public MAdOperatorBase
  {
  public:

    regionRemoveOp(pMesh, DiscreteSF *);
    regionRemoveOp(const regionRemoveOp &);
    ~regionRemoveOp() {}

    operationType type() const { return MAd_RREMOVE; }
  
    void setRegion(pRegion);

    void getCavity(pPList *) const;

    void apply();

  private:

    bool checkConstraints() const;
    bool checkGeometry();
    bool evaluateShapes();
    void evaluateLengths() const;

  private:

    pRegion region;

    bool classifyFaces;
    pGEntity geoFace;

    bool delFace[4];
  };

  // -------------------------------------------------------------------
  inline regionRemoveOp::regionRemoveOp(const regionRemoveOp & _rr):
    MAdOperatorBase(_rr), region(_rr.region),
    classifyFaces(_rr.classifyFaces), geoFace(_rr.geoFace)
  {
    for (int i=0; i<4; i++) delFace[i] = _rr.delFace[i];
  }

  // -------------------------------------------------------------------
  inline regionRemoveOp::regionRemoveOp(pMesh _m, DiscreteSF * _sf):
    MAdOperatorBase(_m,_sf), region(NULL),
    classifyFaces(false), geoFace(NULL)
  {
    for (int i=0; i<4; i++) delFace[i] = false;
  }

  // -------------------------------------------------------------------
  inline void regionRemoveOp::setRegion(pRegion _region)
  {
    region = _region;
  }

  // -------------------------------------------------------------------
  inline void regionRemoveOp::getCavity(pPList * cavity) const
  {
    PList_append(*cavity,(pEntity)region);
  }

  // -------------------------------------------------------------------

}

#endif
