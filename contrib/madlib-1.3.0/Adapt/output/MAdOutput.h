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

#ifndef _H_MADOUTPUT
#define _H_MADOUTPUT

#include "SizeFieldBase.h"
#include "MeshDataBaseInterface.h"

#include <string>
#include <set>

namespace MAd {

  // -------------------------------------------------------------------
  typedef enum MAdOutputData {
    OD_CONSTANT = 0,
    OD_MEANRATIO = 1,
    OD_ORIENTEDMEANRATIO = 2,
    OD_SIZEFIELD_MEAN = 3,
    OD_SIZEFIELD_MIN = 4,
    OD_SIZEFIELD_MAX = 5,
    OD_DIMENSION = 6,
    OD_ITERORDER = 7,
    OD_CURVATURE_DIV = 8,
    OD_CURVATURE_MAX = 9,
    OD_CURVATURE_MIN = 10,
    OD_CURVATURE_MAX_VEC = 11,
    OD_CURVATURE_MIN_VEC = 12,
    OD_ANISO_SF_AXIS0 = 13,
    OD_ANISO_SF_AXIS1 = 14,
    OD_ANISO_SF_AXIS2 = 15,
    OD_DISTANCE       = 16
  } MAdOutputData;

  // -------------------------------------------------------------------

  void MAdGmshOutput(const pMesh m, const pSField sf, 
                     const char *fn, MAdOutputData type);

  void MAdAttachedNodalDataOutput   (const pMesh m, const char *fn, 
                                     pMeshDataId id);
  void MAdAttachedNodalDataVecOutput(const pMesh m, const char *fn, 
                                     pMeshDataId id);

  void printPosEntities(const pPList ents, std::string fn, MAdOutputData type, 
                        const pSField sf=NULL, int id=0);

  void printPosRegions(const std::set<pRegion>, std::string fn, MAdOutputData type, 
                       const pSField sf=NULL, int id=0);

  // -------------------------------------------------------------------

}

#endif
