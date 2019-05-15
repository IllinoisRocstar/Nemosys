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

#ifndef _H_ELEMENTSTATISTICS
#define _H_ELEMENTSTATISTICS

#include "MAdDefines.h"

namespace MAd {

  // -------------------------------------------------------------------
  class ElementStatistics {

  public:

    ElementStatistics();
    ElementStatistics(const ElementStatistics &);
    ~ElementStatistics() {};

    // interface to set information
    void reset();
    void setWorstShape(double v) { worstShape = v; }
    void setMaxLenSq(double v)   { maxLenSq   = v; }
    void setMinLenSq(double v)   { minLenSq   = v; }

    // interface to get information
    double getWorstShape() const { return worstShape; }
    double getMaxLenSq()   const { return maxLenSq;   }
    double getMinLenSq()   const { return minLenSq;   }

  private:

    double worstShape;
    double minLenSq, maxLenSq;
  };

  // -------------------------------------------------------------------
  inline ElementStatistics::ElementStatistics()
  {
    reset();
  }

  // -------------------------------------------------------------------
  inline ElementStatistics::ElementStatistics(const ElementStatistics & eq)
  {
    worstShape   = eq.worstShape;
    minLenSq     = eq.minLenSq;
    maxLenSq     = eq.maxLenSq;
  }

  // -------------------------------------------------------------------
  inline void ElementStatistics::reset()
  {
    worstShape = MAdBIG;
    minLenSq   = MAdBIG;
    maxLenSq   = 0.0;
  }

  // -------------------------------------------------------------------

}

#endif
