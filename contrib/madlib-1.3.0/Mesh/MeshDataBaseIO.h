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
// Authors: Jean-Francois Remacle, Gaetan Compere, Cecile Dobrzynski
// -------------------------------------------------------------------

#ifndef _LOADGMSHMESH_H_
#define _LOADGMSHMESH_H_

#include "MeshDataBaseInterface.h"
#include "MeshDataBaseCommPeriodic.h"

namespace MAd {

  // -------------------------------------------------------------------

  // Load a Gmsh mesh
  //   - msh1 or msh2
  //   - serial or parallel (format msh2 for parallel)
  //   - periodic or non-periodic
  void LoadGmshMesh (MAd::pMesh, const char *);

  // Save a mesh
  //   - msh1 or msh2
  //   - serial only
  //   - If a partitioning table is submitted, 
  //     write the right partition numbers in the file
  void SaveGmshMesh (const MAd::pMesh, const char *, int version=2, 
                     bool saveAll=true, const int * partitionTable=NULL);

  void SaveGmshMeshPer (const MAd::pMesh, const char *,  MDB_DataExchangerPeriodic &deperiodic,int version=1);

  // -------------------------------------------------------------------

}

#endif
