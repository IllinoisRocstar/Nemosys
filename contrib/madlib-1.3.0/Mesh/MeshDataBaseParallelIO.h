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
// Authors: Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef _LOADPARALLELGMSHMESH_H_
#define _LOADPARALLELGMSHMESH_H_

#include "MeshDataBaseInterface.h"

namespace MAd {

  // -------------------------------------------------------------------

#ifdef PARALLEL
  // Save a mesh (for parallel only) \ingroup parallel
  //   - msh1 or msh2
  // GCTODO: merge serial and parallel versions
  void SaveGmshMeshParallel (const pMesh, const char *filename, int version=2);
#endif

  // Load a Gmsh mesh written in msh1 (deprecated)
  // void LoadGmshParallelOld (pMesh, const char *, const int numproc, int version=1);

  // -------------------------------------------------------------------

}

#endif
