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

#ifndef _H_MADFIELDEVALUATORBASE
#define _H_MADFIELDEVALUATORBASE

#include <ostream>
#include <string>

namespace MAd {

  // ----------------------------------------------------------------------
  class MAdFieldEvaluator{
  public:
    virtual bool eval (const double space[3], double time, double* val) const = 0;
    virtual int order () const = 0;
    virtual int nbVar () const = 0;
    virtual ~MAdFieldEvaluator(){};

    virtual void describe(std::ostream& out, const std::string& prefix) const {out << prefix << "Dummy field evaluator \n";} 
  };

  // ----------------------------------------------------------------------
  class MAdConstantScalar : public MAdFieldEvaluator
  {
    double X;
  public :
    MAdConstantScalar(double v):X(v){}
    inline bool eval (const double space[3], double time, double* val) const {
      val[0] = X;
      return true;
    }
    inline int order() const{
      return 0;
    }
    virtual int nbVar () const {return 1;}

    virtual void describe(std::ostream& out,const std::string& prefix="") const{
      out << prefix << "Constant scalar field evaluator u(x,y,z,t) = " << X << "\n";
    }
  
  };

  // ----------------------------------------------------------------------
  class MAdConstantVector : public MAdFieldEvaluator
  {
    double U,V,W;
  public:
    MAdConstantVector(double ux=0,double uy=0,double uz=0):U(ux),V(uy),W(uz){}
    inline bool eval (const double space[3], double time, double* val) const{
      val[0] = U;
      val[1] = V;
      val[2] = W;
      return true;
    }
    inline int order() const{
      return 0;
    }
    virtual int nbVar () const {return 3;}
  
    virtual void describe(std::ostream& out,const std::string& prefix="") const{
      out << prefix << "Constant vector field evaluator u(x,y,z,t) = (" << U << "," << V << "," << W << ")\n";
    }
  
  };

  // ----------------------------------------------------------------------

}

#endif
