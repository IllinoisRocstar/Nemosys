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

#ifndef _H_MEANRATIOEVALUATOR
#define _H_MEANRATIOEVALUATOR

#include "ElementEvaluatorBase.h"

/* -------------------------------------------------------------------

Computes the cubic of the mean ratio of any element of the mesh (tri or tet)
The mean ratio is expressed as:

For tetrahedrons:
     
eta = 12 * ( ( 3*Volume ) ^ (2/3) ) / ( sum(edgeLength^2) )

or equivalently:

eta ^ 3 = 15552 * (Volume) ^ 2 / ( sum(edgeLength^2) ) ^ 3

For triangles:

eta ^ 2 = 48 * ( Area ^ 2 ) / ( sum(edgeLength^2) ) ^ 2

This ensures ratios ranging from 0 (flat element) to 1 (equilateral)

------------------------------------------------------------------- */

namespace MAd {

  // -------------------------------------------------------------------
  class meanRatioEvaluator : public elementEvaluatorBase
  {
  public: 
    meanRatioEvaluator(const DiscreteSF * f);
    virtual ~meanRatioEvaluator() {}

    virtual int F_shape(const pFace pf, const double normal[3], double * result) const;
    virtual int F_shapeWithDisp(const pFace pf, const double normal[3],
                                double disp[3][3], double * result) const;
    virtual int R_shape(const pRegion pr, double * result) const;
    virtual int R_shapeWithDisp(const pRegion rgn, double disp[4][3], 
                                double * shape) const;
    virtual int XYZ_F_shape(const double[3][3], const pMSize[3], const double normal[3], double * result) const;
    virtual int XYZ_R_shape(const double[4][3], const pMSize[4], double * result) const;

    virtual evaluationType getType () const {return MEANRATIO;}
    virtual std::string getName () const {return "Mean ratio";}

  protected:
    virtual int XYZ_F_shape(const double[3][3], const pMSize pS, const double normal[3], double * result) const;
    virtual int XYZ_R_shape(const double[4][3], const pMSize pS, double * result) const;
  };

  // -------------------------------------------------------------------

}

#endif
