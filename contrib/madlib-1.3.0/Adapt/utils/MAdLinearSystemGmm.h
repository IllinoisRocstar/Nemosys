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

#ifndef _H_LINEARSYSTEMGMM_MAD
#define _H_LINEARSYSTEMGMM_MAD

#ifdef _HAVE_GMM_

// Interface to GMM++

#include "MAdLinearSystem.h"

#include <gmm.h>

namespace MAd {

  // -------------------------------------------------------------------
  class MAdLinearSystemGmm : public MAdLinearSystem {

  public :

    MAdLinearSystemGmm () : MAdLinearSystem(), _a(0), _b(0), _x(0), _prec(1.e-8), _noisy(0) {}
    ~MAdLinearSystemGmm ()
    {
      delete _a;
      delete _b;
      delete _x;
    }

  public:

    bool isAllocated () const {return _a != 0;}
    void allocate (int _nbRows)
    {
      if (_a) delete _a;
      if (_x) delete _x;
      if (_b) delete _b;
      _a = new  gmm::row_matrix< gmm::wsvector<double> >(_nbRows,_nbRows);
      _b = new  std::vector<double>(_nbRows);
      _x = new  std::vector<double>(_nbRows);    
    }
    void  addToMatrix    (int _row, int _col, double _val) 
    {
      if (_val != 0.0) (*_a)(_row, _col) += _val;
    }
    double getFromMatrix (int _row, int _col) const
    {
      return (*_a)(_row, _col);
    }
    void  addToRightHandSide    (int _row, double _val) 
    {
      if (_val != 0.0) (*_b)[_row]+=_val;
    }
    double getFromRightHandSide (int _row) const 
    {
      return (*_b)[_row];
    }
    double getFromSolution (int _row) const 
    {
      return (*_x)[_row];
    }
    void zeroMatrix () 
    {
      gmm::clear(*_a);
    }
    void zeroRightHandSide () 
    {
      for (unsigned int i = 0; i < _b->size(); i++) (*_b)[i] = 0;
    }
    void reorder () {
    
    }
    void setPrec(double p){_prec=p;}
    void setNoisy(int n){_noisy=n;}
    int systemSolve () 
    {
      //gmm::ilut_precond< gmm::row_matrix< gmm::wsvector<double> > > P(*_a, fill, _prec);
      gmm::ildltt_precond< gmm::row_matrix< gmm::wsvector<double> > > P(*_a, fill, _prec);
      gmm::iteration iter(eps);
      iter.set_noisy(_noisy);
      //     gmm::sequential_additive_schwarz(*_a, *_x, *_b, P, vB, iter, local_solver, global_solver)
      if (sType==GMRES)gmm::gmres(*_a, *_x, *_b, P, 100, iter);  // execute the GMRES algorithm
      else gmm::cg(*_a, *_x, *_b, P, iter);  // execute the CG algorithm
      return 1;
    }

  private:

    gmm::row_matrix<gmm::wsvector<double> > *_a;
    std::vector<double> *_b, *_x;
    double _prec;
    int _noisy;
  };

  // -------------------------------------------------------------------

#endif

}

#endif
