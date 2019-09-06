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
// Authors: Cecile Dobrzynski, Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef _MARK_H_
#define _MARK_H_

#include "MeshDataBaseInterface.h"

#include <vector>

namespace MAd {

  // -------------------------------------------------------------------
#ifdef PARALLEL
  void MarkTriangles       (pMesh mesh, pMeshDataId tagElt);
  void MarkTrianglesSmooth (pMesh mesh, pMeshDataId tagElt);
  void MarkTrianglesRandom (pMesh mesh, pMeshDataId tagElt);
  void MarkTets            (pMesh mesh, pMeshDataId tagElt);
  void MarkTetsSmooth      (pMesh mesh, pMeshDataId tagElt);
  void MarkTetsRandom      (pMesh mesh, pMeshDataId tagElt);
  int  MarkTetsManifold    (pMesh mesh, pMeshDataId tagElt);
#endif

  // -------------------------------------------------------------------
  void VertexTagMove         (pMesh mesh, pMeshDataId tagMove);
  void VertexTagMoveList     (pMesh mesh, std::vector<std::vector<int> >& transfo,
                              pMeshDataId tagMove, pMeshDataId tagTransfo);
  void VertexTagMoveInverse  (pMesh mesh, pMeshDataId tagMove);

  // -------------------------------------------------------------------
  void MarkPeriodicTriangles (pMesh mesh, std::vector<std::vector<int> >& transfo,
                              pMeshDataId tagElt, pMeshDataId tagMove, pMeshDataId tagTransfo);
  void MarkPeriodicTets      (pMesh mesh, std::vector<std::vector<int> >& transfo,
                              pMeshDataId tagElt, pMeshDataId tagMove, pMeshDataId tagTransfo);
  int  MarkGroupPeriodicTets (pMesh mesh, pMeshDataId tagElt, pMeshDataId tagMove);

  // -------------------------------------------------------------------

}

#endif
