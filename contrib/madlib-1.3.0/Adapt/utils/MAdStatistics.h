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

#ifndef _H_MADSTATISTICS
#define _H_MADSTATISTICS

#include "MAdSingleton.h"

#include <iostream>

namespace MAd {

  // -------------------------------------------------------------------
  class MAdStatistics {

  public:
  
    MAdStatistics()  {};
    ~MAdStatistics() {};
  
    void initialize();
    void finalize();

    void print(std::ostream& out) const;

    // -----------------------
    // CPU time
    // -----------------------

  private:

    double t_eSplits, t_eCollapses, t_eSwaps, t_rSlivers, t_fSlivers;

  public:
  
    void addCPUESplits    (double dt) { t_eSplits    += dt; }
    void addCPUECollapses (double dt) { t_eCollapses += dt; }
    void addCPUESwaps     (double dt) { t_eSwaps     += dt; }
    void addCPURSlivers   (double dt) { t_rSlivers   += dt; }
    void addCPUFSlivers   (double dt) { t_fSlivers   += dt; }

    double getCPUESplits    () const { return t_eSplits;    }
    double getCPUECollapses () const { return t_eCollapses; }
    double getCPUESwaps     () const { return t_eSwaps;     }
    double getCPUSlivers    () const { return ( t_rSlivers + t_fSlivers ); }


    // -----------------------------
    // Elementary operations count
    // -----------------------------

  private:

    int num_eSplits, num_eCollapses, num_eSwaps;

  public:

    void addNumESplits    (int num) { num_eSplits    += num; }
    void addNumECollapses (int num) { num_eCollapses += num; }
    void addNumESwaps     (int num) { num_eSwaps     += num; }
    int getNumESplits     () { return num_eSplits; }
    int getNumECollapses  () { return num_eCollapses; }
    int getNumESwaps      () { return num_eSwaps; }


    // -----------------------------
    // Others
    // -----------------------------

  private:

    int numInfLoops;

  public:

    void addInfiniteLoops(int num) { numInfLoops += num; }

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<MAdStatistics> MAdStatisticsSgl;

}

#endif
