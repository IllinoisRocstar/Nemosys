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

#ifndef _H_MADTIMEMANAGER
#define _H_MADTIMEMANAGER

#include "MAdSingleton.h"

namespace MAd {

  // -------------------------------------------------------------------
  class MAdTimeManager {

  public:
  
    MAdTimeManager() {};
    ~MAdTimeManager() {};
  
    void   initialize();
    void   finalize();

    void   setTime(double);
    void   incrementTime(double);

    double getTime() const;
  
  private:

    double time;

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<MAdTimeManager> MAdTimeManagerSgl;

}

#endif
