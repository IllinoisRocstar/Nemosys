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
  Store a list of pointers to geometric entities

  The list is stored as a standart vector:
    - allow fast random access
    - allow fast iterating
    - can be allocated to a particular size
    - slow at randomly adding/removing elements but not useful here
*/

#ifndef _H_PGLIST
#define _H_PGLIST

#include "ModelInterface.h"

#include <vector>

namespace MAd {

  // -------------------------------------------------------------------
  class PGList {

  public:

    PGList() {};
    PGList(const PGList& ori)
    {
      entities = ori.entities;
    };
  
    ~PGList() {};

  public:

    void clear() { entities.clear(); };

  public:

    std::vector<MAdGEntity *> entities;

  };

  // -------------------------------------------------------------------

}

#endif
