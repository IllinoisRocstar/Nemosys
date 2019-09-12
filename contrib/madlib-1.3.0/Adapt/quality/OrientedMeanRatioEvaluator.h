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

#ifndef _H_ORIENTEDMEANRATIOEVALUATOR
#define _H_ORIENTEDMEANRATIOEVALUATOR

#include "MeanRatioEvaluator.h"

/* -------------------------------------------------------------------

Computes the cubic of the mean ratio of any element of the 
mesh (tri or tet) penalised (for every edge 'e') by a factor 

   | |cos(alpha_e)| - 0.5 | + 0.5

where 'alpha_e' is the angle of 'e' with the direction of the minimal 
size. This is equivalent to the mean ratio for isotropic sizes.

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
  class orientedMeanRatioEvaluator : public meanRatioEvaluator
  {
  public: 
    orientedMeanRatioEvaluator(const DiscreteSF * f);
    virtual ~orientedMeanRatioEvaluator() {}

    virtual evaluationType getType () const {return ORIENTEDMEANRATIO;}
    virtual std::string getName () const {return "Oriented mean ratio";}

  protected:
    virtual int XYZ_F_shape(const double[3][3], const pMSize pS, const double normal[3], double * result) const;
    virtual int XYZ_R_shape(const double[4][3], const pMSize pS, double * result) const;
  };

  // -------------------------------------------------------------------

}

#endif
