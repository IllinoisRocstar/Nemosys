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

#ifndef _PARALLELUTILS_H_
#define _PARALLELUTILS_H_

namespace MAd {

#ifdef _HAVE_METIS_
  void PartitionMesh(pMesh m, int nbPart, const char *);
#endif

#ifdef PARALLEL
  void E_facesAttachId(pMesh mesh,pMeshDataId tagAdj);
  void F_regionsAttachId(pMesh mesh,pMeshDataId tagAdj);
#endif

}

#endif
