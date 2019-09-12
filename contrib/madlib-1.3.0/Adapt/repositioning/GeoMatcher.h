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

#ifndef _H_GEOMATCHER
#define _H_GEOMATCHER

#include "MAdElasticityOp.h"
#include "SliverFaceHandler.h"
#include "SliverRegionHandler.h"
#include <list>

namespace MAd {

  class MeshAdapter;

  // -------------------------------------------------------------------
  class geoMatcher: public MAdElasticityOp
  {

  public:

    geoMatcher(pMesh, MeshAdapter *);
    geoMatcher(const geoMatcher&);
    ~geoMatcher();

    void setForceRelocation(bool _force) { force = _force; }
    void setStrictChecking (bool sc) { strictChecking = sc; }

    void setSliverFaceHandler(sliverFaceHandler * s) { sliverFOp = s; }
    void setSliverRegionHandler(sliverRegionHandler * s) { sliverROp = s; }

    bool snap();

    void reportFailure(double ratio) { failures.push_back(ratio); }
    void printFailures(std::ostream&) const;

  private:

    bool force; // if the nodes are snapped without any check nor repositioning inside the domain
    bool strictChecking; // true: snap() returns true only if we the full snapping is done

    std::list<double> failures;

    MeshAdapter * adapter;

    sliverFaceHandler    *  sliverFOp;
    sliverRegionHandler  *  sliverROp;

  };

  // -------------------------------------------------------------------

}

#endif
