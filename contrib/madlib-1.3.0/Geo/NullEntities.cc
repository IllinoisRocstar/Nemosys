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


#ifndef _HAVE_GMSH_

#include "NullModel.h"

using namespace MAd;

// -------------------------------------------------------------------
void NullGEntity::setPhysical(int d, int t)
{
  //       if ( phys ) phys.delGE(this);
  phys = model->getPhysical(d,t);
  //       phys.addGE(this);
}

// -------------------------------------------------------------------

#endif
