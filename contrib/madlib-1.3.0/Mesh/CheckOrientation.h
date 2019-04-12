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

#ifndef _CHECKORIENTATION_H_
#define _CHECKORIENTATION_H_

#ifdef PARALLEL

#include "MeshDataBaseInterface.h"

namespace  MAd {
  int CheckEdgesOrientation(pMesh mesh);
  int CheckFacesOrientation(pMesh mesh);
}

#endif

#endif
