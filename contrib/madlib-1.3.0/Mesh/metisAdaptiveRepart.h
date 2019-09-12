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

#ifdef PARALLEL
#ifndef H_METISADAPTIVEREPART
#define H_METISADAPTIVEREPART

#include "MeshDataBaseInterface.h"

namespace MAd {

#ifdef _HAVE_PARMETIS_
  //! Repartition (adapt) the graph of the elements over the 
  //! set of processors and attach the destination to the 
  //! elements under the tag 'tagElt'  \ingroup parallel
  extern void metisAdaptiveRepart(pMesh mesh, pMeshDataId tagElt);
#endif

}

#endif
#endif
