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

#ifndef _H_PHYSICAL
#define _H_PHYSICAL

#include <set>

// #ifdef _HAVE_GMSH_
// #include "GmshEntities.h"
// #else
// #include "NullEntities.h"
// #endif

namespace MAd {

  // -------------------------------------------------------------------
  class Physical {
    
  private:
    int _tag;
    int _dim;
//     std::set<MAdGEntity *> gents;

  public:
    
    Physical(int d, int t) { _dim = d; _tag = t; }
    ~Physical() {}
    
//     void addGE(MAdGEntity * newg) { gents.insert(newg); }
//     void delGE(MAdGEntity * delg)
//     {
//       std::set<MAdGEntity *>::iterator it = gents.find(delg);
//       if ( it != gents.end() ) gents.erase(it);
//     }

    int tag() const { return _tag; }
    int dim() const { return _dim; }
    
  };

  // -------------------------------------------------------------------
}

#endif

