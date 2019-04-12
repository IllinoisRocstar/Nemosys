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

#include "MAdModel.h"

namespace MAd {

  Physical * MAdModel::getPhysical(int d, int t)
  {
    std::set<Physical *>::const_iterator it = physicals.begin();
    for (; it != physicals.end(); it++)
      {
        if ( (*it)->dim() == d && (*it)->tag() == t ) return *it;
      }

    Physical * newp = new Physical(d,t);
    physicals.insert(newp);
    return newp;
  }

}
