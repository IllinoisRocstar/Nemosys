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

#ifndef _H_LINEARSYSTEM_MAD
#define _H_LINEARSYSTEM_MAD

namespace MAd {

  // A class that encapsulates a linear system solver interface :
  // building a sparse matrix, solving a linear system

  // -------------------------------------------------------------------
  enum SolverType { GMRES, CG };

  // -------------------------------------------------------------------
  class MAdLinearSystem {

  public:

    MAdLinearSystem (): sType(GMRES),fill(0),eps(1.e-12) {}
    virtual ~MAdLinearSystem () {}

  public:

    virtual bool isAllocated () const = 0;
    virtual void allocate (int nbRows) = 0;
    virtual void  addToMatrix    (int _row, int _col, double val) = 0;
    virtual double getFromMatrix (int _row, int _col) const = 0;
    virtual void  addToRightHandSide    (int _row, double val) = 0;
    virtual double getFromRightHandSide (int _row) const = 0;
    virtual double getFromSolution (int _row) const = 0;
    virtual void zeroMatrix () = 0;
    virtual void zeroRightHandSide () = 0;
    virtual void reorder() = 0;
    virtual int systemSolve () = 0;
		virtual void set_nnz(int row,int nz){}
		virtual void allocate_matrix(){}

    virtual void setSolver(SolverType _type) {sType = _type;}
    void setFillIn(int _k) {fill=_k;}
    void setEps(double _eps) {eps=_eps;}
    virtual void setPrec(double){}
    virtual void setNoisy(int){}

  protected:

    SolverType sType;
    int fill;
    double eps;

  };

  // -------------------------------------------------------------------
  class MAdLinearSystemNull : public MAdLinearSystem {

  public:
  
    MAdLinearSystemNull ()
    {
      printf("Error: no linear systems solver available\n");
      throw;
    }

  public:

    bool isAllocated () const { return false; }
    void allocate (int nbRows) {}
    void  addToMatrix    (int _row, int _col, double val) {}
    double getFromMatrix (int _row, int _col) const { return 0.; }
    void  addToRightHandSide    (int _row, double val) {}
    double getFromRightHandSide (int _row) const { return 0.; }
    double getFromSolution (int _row) const { return 0.; }
    void zeroMatrix () {}
    void zeroRightHandSide () {}
    void reorder() {}
    int systemSolve () { return 0; }
  };

  // -------------------------------------------------------------------

}

// -------------------------------------------------------------------

#ifdef _HAVE_PETSC_
#include "MAdLinearSystemPETSc.h"
namespace MAd {
  typedef MAdLinearSystemPETSc MAdLinearSystemDef;
}
#else

#ifdef _HAVE_SPARSKIT_
#include "MAdLinearSystemSparskit.h"
namespace MAd {
  typedef MAdLinearSystemSparskit MAdLinearSystemDef;
}
#else

#ifdef _HAVE_GMM_
#include "MAdLinearSystemGmm.h"
namespace MAd {
  typedef MAdLinearSystemGmm MAdLinearSystemDef;
}
#else
namespace MAd {
  typedef MAdLinearSystemNull MAdLinearSystemDef;
}
#endif // end gmm
#endif // end sparskit
#endif // end PETSc

// -------------------------------------------------------------------

#endif


