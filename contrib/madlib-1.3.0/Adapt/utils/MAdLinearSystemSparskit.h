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

#ifndef _H_LINEARSYSTEMSPARSKIT_MAD
#define _H_LINEARSYSTEMSPARSKIT_MAD

// -------------------------------------------------------------------
//  Interface to Y Saad's sparskit lib
// -------------------------------------------------------------------

#ifdef _HAVE_SPARSKIT_

#include "MAdLinearSystem.h"
#include "MAdMessage.h"

#include "CSR_Matrix.h"
#include "CSR_Vector.h"

namespace MAd {

  // -------------------------------------------------------------------
  class MAdLinearSystemSparskit : public MAdLinearSystem {

  public:

    MAdLinearSystemSparskit ():
      MAdLinearSystem(),_a(0),_b(0),_x(0)
    {
    }
    ~MAdLinearSystemSparskit ()
    {
      delete _a;
      delete _b;
      delete _x;
    }
    bool isAllocated () const {return _a != 0;}
    void allocate (int _nbRows)
    {
      if (_a) delete _a;
      if (_b) delete _b;
      if (_x) delete _x;
      _a = new CSR_Matrix (_nbRows);
      _b = new CSR_Vector (_nbRows);
      _x = new CSR_Vector (_nbRows);
    }
    void addToMatrix (int _row, int _col, double _val) {
      if (_val != 0.0) _a->AddMatrix (_row+1, _col+1, _val);
    }
    double getFromMatrix (int _row, int _col) const {
      return _a->GetMatrix (_row+1,_col+1);
    }
    void addToRightHandSide (int _row, double _val) {
      if (_val != 0.0) _b->AddVal(_row+1, _val);
    }
    double getFromRightHandSide (int _row) const {
      return _b->GetVal (_row+1);
    }
    double getFromSolution (int _row) const {
      return _x->GetVal (_row+1);
    }
    void zeroMatrix () {
      _a->ZeroMatrix();
    }
    void zeroRightHandSide () {
      _b->ZeroArray();
    }
    void reorder () {
    
    }
    void setSolver(SolverType _type) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"CG not available in Sparskit, GMRES will be used\n");
      sType = GMRES;
    }
    int systemSolve () {
      _a->EndOfAssembly();
      std::string solver = "gmres";
      //    if (sType == CG) solver = "cg";
      SPARSKIT_LINEAR_SOLVER_ ( "rcmk","ilut",solver,fill,eps,*_a,*_b,*_x);
      return 1;
    }

  private:
  
    CSR_Matrix *_a;
    CSR_Vector *_b;
    CSR_Vector *_x;
  
  };

  // -------------------------------------------------------------------

}

#endif

#endif
