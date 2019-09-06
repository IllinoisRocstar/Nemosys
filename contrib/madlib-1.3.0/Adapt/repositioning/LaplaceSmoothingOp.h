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

#ifndef _H_LAPLACESMOOTHINGOP
#define _H_LAPLACESMOOTHINGOP

#include "NodesRepositioningOp.h"

namespace MAd {

  // -------------------------------------------------------------------
  class LaplaceSmoothingOp : public nodesRepositioningOp {

  public:

    LaplaceSmoothingOp(pMesh m, DiscreteSF * sf);
    ~LaplaceSmoothingOp() {}

  public:

    virtual int run(double * L2Norm);
    virtual int runOptimal(double * L2Norm);
    virtual int runFast(double * L2Norm);

  };

  // -------------------------------------------------------------------

}

#endif
