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

#ifndef _H_OPERATORTOOLS
#define _H_OPERATORTOOLS

#include "MeshDataBaseInterface.h"

namespace MAd {

  // -------------------------------------------------------------------

  bool V_setPosition(pVertex, double[3]);

  // --- Related to edge collapse ---

  void E_collapse(pMesh,   // the mesh
                  pEdge,   // the edge to be collapsed
                  pVertex, // the vertex to be deleted
                  pVertex  // the target vertex
                  );

  void E_collapseOnGFace(pMesh,   // the mesh
                         pEdge,   // the edge to be collapsed
                         pVertex, // the vertex to be deleted
                         pVertex  // the target vertex
                         );

  // --- Related to edge split ---

  pVertex E_split(pMesh, pEdge, double[3], double t=0.5); 

  // -------------------------------------------------------------------

}

#endif

