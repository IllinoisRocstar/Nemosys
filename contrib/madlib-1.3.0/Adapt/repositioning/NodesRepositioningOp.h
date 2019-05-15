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

#ifndef _H_NODESREPOSITIONINGOP
#define _H_NODESREPOSITIONINGOP

#include "DiscreteSF.h"
#include "MeshQualityManager.h"
#include "MeshDataBaseInterface.h"

namespace MAd {

  // -------------------------------------------------------------------
  enum repositionType {
    LAPLACE_FAST,
    LAPLACE_OPTIMAL
  };

  // -------------------------------------------------------------------
  class nodesRepositioningOp {

  public:

    nodesRepositioningOp(pMesh m, DiscreteSF * sf): 
      mesh(m), sizeField(sf), mqm(MeshQualityManagerSgl::instance())
    {
      dim = 3;
      if ( M_numRegions(mesh) == 0 ) dim = 2;
    }
    virtual ~nodesRepositioningOp() {}

  public:

    virtual int run(double * L2Norm)=0;

  protected:

    pMesh mesh;
    int dim;
    DiscreteSF * sizeField;
    MeshQualityManager& mqm;

  };

  // -------------------------------------------------------------------

}

#endif
