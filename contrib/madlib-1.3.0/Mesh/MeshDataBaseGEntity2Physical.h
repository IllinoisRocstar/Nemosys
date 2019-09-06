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

#ifndef H_MESHDATABASEGENTITY2PHYSICAL
#define H_MESHDATABASEGENTITY2PHYSICAL
#include <map>
#include "ModelInterface.h"
// build a reverse map of geomFeatures_Tags
// (a map would be enough as we read only one physical tag for now,
// but a multimap is used in prevision...)

namespace MAd {

  class GEntity2Physical:private std::multimap<pGEntity,int>{
  public:
    GEntity2Physical(std::multimap<int,pGEntity> &inmap);
    int get_first_tag(pGEntity g);
  };
  
}

#endif
