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

#include "LaplaceSmoothingOp.h"
#include "VertexMoveOp.h"
#include "MeshParametersManager.h"
#include "MathUtils.h"

#include <vector>
using std::vector;
#include <math.h>

namespace MAd {

  // -------------------------------------------------------------------
  LaplaceSmoothingOp::LaplaceSmoothingOp(pMesh m, DiscreteSF * sf): 
    nodesRepositioningOp(m,sf)
  {}

  // -------------------------------------------------------------------
  int LaplaceSmoothingOp::run(double * L2Norm)
  {
    return runOptimal(L2Norm);
  }

  // -------------------------------------------------------------------
  int LaplaceSmoothingOp::runOptimal(double * L2Norm)
  {
    *L2Norm = 0.0;

    vertexMoveOp * vMove = new vertexMoveOp(mesh,sizeField,true);

    VIter vit = M_vertexIter(mesh);
    while( pVertex pv = VIter_next(vit) ) {

      if ( V_whatInType(pv) != dim ) continue;

      double ori[3];
      V_coord(pv,ori);

      double opt[3];
      if (vMove->computeOptimalLocation(pv, opt)) vMove->setPosition(pv, opt);

      double worst;
      if( vMove->evaluate(&worst) ) {
        vMove->apply();
        double disp[3];
        diffVec(ori,opt,disp);
        *L2Norm += dotProd(disp,disp);
      }

      vMove->resetDisplacements();
    }

    return 1;
  }

  // -------------------------------------------------------------------
  int LaplaceSmoothingOp::runFast(double * L2Norm)
  {
    *L2Norm = 0.0;

    vertexMoveOp * vMove = new vertexMoveOp(mesh,sizeField,true);

    VIter vit = M_vertexIter(mesh);
    while( pVertex pv = VIter_next(vit) ) {

      if ( V_whatInType(pv) != dim ) continue;

      double ori[3];
      V_coord(pv,ori);

      double opt[3];
      V_cavityCenter(pv,opt);
      vMove->setPosition(pv, opt);

      double worst;
      if( vMove->evaluate(&worst) ) {
        vMove->apply();
        double disp[3];
        diffVec(ori,opt,disp);
        *L2Norm += dotProd(disp,disp);
      }

      vMove->resetDisplacements();
    }

    return 1;
  }

  // -------------------------------------------------------------------

}

