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

#include "MeshDataBaseGEntity2Physical.h"

namespace MAd {

GEntity2Physical::GEntity2Physical(std::multimap<int,pGEntity> &inmap){
  for(std::multimap<int,pGEntity>::iterator it=inmap.begin(); it!=inmap.end(); it++){
    insert(std::pair<pGEntity,int>(it->second,it->first));
  }
}

int GEntity2Physical::get_first_tag(const pGEntity g){
  iterator first_it,end_it;
  first_it=lower_bound(g);
  end_it=upper_bound(g);
  if(first_it==end_it){
    // This entity doesn't have physical tag
    // I don't think this could happen but who knows...
    return 0;
  }
  return first_it->second;
}

}

