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

/*
  Store a list of pointers to entities

  The list is stored as a standart vector:
    - allow fast random access
    - allow fast iterating
    - can be allocated to a particular size
    - slow at randomly adding/removing elements but not useful here
*/

#ifndef _H_PLIST
#define _H_PLIST

#include "MeshDataBase.h"

#include <vector>

namespace MAd {

  // -------------------------------------------------------------------
  class PList {

  public:

    PList() {};
    PList(const PList& ori)
    {
      entities = ori.entities;
    };
  
    ~PList() {};

  public:

    void clear() { entities.clear(); };

  public:

    std::vector<MDB_MeshEntity *> entities;

  };

  // -------------------------------------------------------------------

}
#endif
