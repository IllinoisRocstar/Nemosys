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

#ifndef _H_CONSTRAINT
#define _H_CONSTRAINT

#include "MeshDataBaseInterface.h"

namespace MAd {

  // -------------------------------------------------------------------

  // --- entity level constraining ---

  void EN_constrain   (pEntity);
  void EN_unconstrain (pEntity);
  bool EN_constrained (pEntity);

  // --- mesh level constraining ---

#ifdef PARALLEL
  // Constrain parallel interface elements 
  void UpdateParallelConstraint(pMesh mesh);
  // unConstrain parallel interface elements
  void DeleteParallelConstraint(pMesh mesh);
#endif

  // Constrain periodic interface elements
  void UpdatePeriodicConstraint3d(pMesh mesh);
  void UpdatePeriodicConstraint2d(pMesh mesh);

  // Remove all constraints on elements (periodic, parallel and misc.)
  void DeleteConstraint(pMesh mesh);

  // -------------------------------------------------------------------

}

#endif
