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

#ifndef _H_CALLBACKDEFINITION
#define _H_CALLBACKDEFINITION

#include "MeshDataBaseInterface.h"
#include "MAdOperations.h"

namespace MAd {

  // -------------------------------------------------------------------
  // callback functions
  typedef void (*CBFunction)(pPList, pPList, 
                             void *, operationType, pEntity);
  typedef void (*CBFunc_move)(pVertex, double *, void *);

  // -------------------------------------------------------------------
  struct CBStructure {

    CBFunction function;
    void* data;

    bool operator== (CBStructure cb) const {
      if ( function==cb.function && data==cb.data) return true;
      return false;
    }
  };

  // -------------------------------------------------------------------
  struct CBStructureMove {
    CBFunc_move function;
    void* data;

    bool operator== (CBStructureMove cbm) const {
      if ( function==cbm.function && data==cbm.data) return true;
      return false;
    }
  };

}

#endif
