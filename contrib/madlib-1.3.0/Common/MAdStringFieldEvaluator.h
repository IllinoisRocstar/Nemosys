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

#ifndef _H_MADSTRINGFIELDEVALUATOR
#define _H_MADSTRINGFIELDEVALUATOR

#include "MAdFieldEvaluator.h"

#ifdef _HAVE_MATHEX_

#include "mathex.h"

#include <list>
#include <cstdarg>
#include <string>
#include <vector>

namespace MAd {

  // ----------------------------------------------------------------------
  class MAdStringFieldEvaluator : public MAdFieldEvaluator
  {

  public:

    ~MAdStringFieldEvaluator ();
    // int numberOfStrings is the number of "const char *" given in "..."
    MAdStringFieldEvaluator (int numberOfStrings, ...);
    MAdStringFieldEvaluator(const std::vector<std::string>&);

  private:

    void buildEvaluators(const std::vector<std::string>& exp);

  public:

    bool eval (const double space[3], double time, double * val) const ;
    bool eval (const std::vector<double> space, double time, double * val) const ;

    int nbVar() const;
    int order() const;

  private:

    mutable double x,y,z,t;
    std::list<smlib::mathex*> expressions;

  };

}

#else

#include "MAdMessage.h"

namespace MAd {

  // ----------------------------------------------------------------------
  class MAdStringFieldEvaluator : public MAdFieldEvaluator
  {

  public:

    ~MAdStringFieldEvaluator () {}
    MAdStringFieldEvaluator (int, ...)
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "No string evaluation library available");
    }
    MAdStringFieldEvaluator(const std::vector<std::string>&)
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "No string evaluation library available");
    }

  private:

    void buildEvaluators(const std::vector<std::string>&) { throw; }

  public:

    bool eval (const double space[3], double time, double * val) const 
    { throw; return false; }
    bool eval (const std::vector<double> space, double time, double * val) const 
    { throw; return false; }

    int nbVar() const { throw; return -1; }
    int order() const { throw; return -1; }

  };

  // ----------------------------------------------------------------------

}

#endif

#endif
