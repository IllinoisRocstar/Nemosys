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

#ifndef _H_MADRESOURCEMANAGER
#define _H_MADRESOURCEMANAGER

#include "MAdSingleton.h"
#include <iostream>
#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  class MAdResourceManager {

  public:
  
    MAdResourceManager()  { reset(); };
    ~MAdResourceManager() {};
  
    void initialize();
    void finalize();

    void reset() { initialTime = getTime(); };
    double elapsedTime() const;
    double getTime() const;

    //! Gives the amount of memory currently used in the RAM (resident_size) 
    //! and in total (RAM, cache, disk) (virtual_size) \ingroup tools
    bool getMemoryUsage(double& resident_size, double& virtual_size) const;
    void printMemoryUsage(std::string step="", std::ostream& out=std::cout) const;

  private:

    double initialTime;

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<MAdResourceManager> MAdResourceManagerSgl;
}

#endif
