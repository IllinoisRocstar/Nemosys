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

#include "MAdTimeManager.h"

namespace MAd {

  // -------------------------------------------------------------------
  void MAdTimeManager::initialize()
  {
    time = 0.;
  }

  // -------------------------------------------------------------------
  void MAdTimeManager::finalize()
  {
  }

  // -------------------------------------------------------------------
  void MAdTimeManager::setTime(double _t)
  {
    time = _t;
  }

  // -------------------------------------------------------------------
  void MAdTimeManager::incrementTime(double dt)
  {
    time += dt;
  }

  // -------------------------------------------------------------------
  double MAdTimeManager::getTime() const
  {
    return time;
  }

  // -------------------------------------------------------------------

}
