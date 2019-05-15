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

#ifdef _HAVE_MATHEX_

#include "MAdStringFieldEvaluator.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
using std::list;
using std::vector;
using std::string;

namespace MAd {

  // ----------------------------------------------------------------------
  MAdStringFieldEvaluator :: MAdStringFieldEvaluator (int numberOfStrings, ...):
    MAdFieldEvaluator() 
  {
    vector<string> exp;
  
    va_list ap;
    va_start(ap,numberOfStrings);
    for (int i=0;i<numberOfStrings;i++){
      const char *s = va_arg(ap,const char *); 
      string str = s;
      exp.push_back(str);
    }
    va_end(ap);

    buildEvaluators(exp);
  }

  // ----------------------------------------------------------------------
  MAdStringFieldEvaluator::MAdStringFieldEvaluator(const vector<string>& exp):
    MAdFieldEvaluator()
  {
    buildEvaluators(exp);
  }

  // ----------------------------------------------------------------------
  void MAdStringFieldEvaluator::buildEvaluators(const vector<string>& exp)
  {
    vector<string>::const_iterator eIter = exp.begin();
    for (;eIter!=exp.end();eIter++) {

      smlib::mathex * expr = new smlib::mathex();
      expr->addvar("x",&x);
      expr->addvar("y",&y);
      expr->addvar("z",&z);
      expr->addvar("t",&t);

      bool success = false;
      bool parsed  = false;
      try {
        expr->expression(*eIter);
        expr->parse();
        parsed  = true;
        success = true;
      }
      catch(smlib::mathex::error e) {
        std::cout << e.what() << std::endl;
        if(!parsed) {
          std::cout << *eIter << std::endl;
          std::cout << std::setw(expr->stopposition()) << "^" << std::endl;
        }
      }
      if ( !success ){
        std::cerr << "Error when parsing string field evaluator: possible symbols are exp, log, log10, sqrt, sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, abs, round, floor, ceil, trunc, deg, fac, rad\n";
        exit(0);
      }
      expressions.push_back(expr);
    }
  }

  // ----------------------------------------------------------------------
  MAdStringFieldEvaluator::~MAdStringFieldEvaluator()
  {
    std::list<smlib::mathex*>::iterator it = expressions.begin();
    for (; it!=expressions.end(); ++it) {
      delete(*it);
    }
  }

  // ----------------------------------------------------------------------
  bool MAdStringFieldEvaluator::eval(const vector<double> space, 
                                            double time, 
                                            double * val) const 
  {
    if (space.size() != 3) throw;

    double spaceTbl[3];
    for (int i=0; i<3; i++) spaceTbl[i] = space[i];

    return eval(spaceTbl, time, val);
  }

  // ----------------------------------------------------------------------
  bool MAdStringFieldEvaluator::eval(const double space[3], 
                                            double time, 
                                            double * val) const 
  {
    x = space[0];
    y = space[1];
    z = space[2];
    t = time;
  
    int i = 0;
    std::list<smlib::mathex*>::const_iterator it = expressions.begin();
    for (; it!=expressions.end(); ++it) {
      val[i++] = (*it)->eval();
    }

    return true;
  }

  // ----------------------------------------------------------------------
  int MAdStringFieldEvaluator::nbVar() const
  {
    return expressions.size();
  }

  // ----------------------------------------------------------------------
  int MAdStringFieldEvaluator::order() const
  {
    return 3;
  }

  // ----------------------------------------------------------------------

}

#endif
