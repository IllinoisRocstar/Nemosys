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

#ifndef _CHECKMESH_H_
#define _CHECKMESH_H_

#include "madlib_export.h"

#include "MeshDataBaseInterface.h"

namespace MAd {

  // -------------------------------------------------------------------
  typedef enum MeshStatus { 
    VALID                   = 0,
    NEGATIVE_VOLUME         = 1,
    GEOM_INCOMPATIBILITY    = 2,
    WRONG_EDGE_TO_RGN_CONN  = 3,
    WRONG_FACE_TO_RGN_CONN  = 4,
    WRONG_FACE_TO_VTX_CONN  = 5,
    WRONG_EDGE_CONN         = 6,
    WRONG_ENTITY_POINTERS   = 7,
    WRONG_ITERATORS         = 8,
    WRONG_PARAMETERS        = 9
  } MeshStatus;

  // -------------------------------------------------------------------
  typedef enum checkType { 
    CHECK_ALL,
    CHECK_VOLUME,
    CHECK_GEOM_COMPATIBILITY,
    CHECK_EDGE_TO_RGN_CONN,
    CHECK_FACE_TO_RGN_CONN,
    CHECK_FACE_TO_VTX_CONN,
    CHECK_EDGE_CONN,
    CHECK_ENTITY_POINTERS,
    CHECK_ITERATORS,
    CHECK_PARAMETERS
  } checkType;

  // -------------------------------------------------------------------
  // return 1 if the mesh passes the test successfully, 0 otherwise
  bool MADLIB_EXPORT checkMesh(MDB_Mesh *mesh,
                               checkType type=CHECK_ALL,
                               int verbose=1,
                               std::ostream &out=std::cout,
                               MeshStatus *status=NULL);

  // -------------------------------------------------------------------

}

#endif
