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

#ifndef _H_ELEMENTEVALUATORBASE
#define _H_ELEMENTEVALUATORBASE

#include "DiscreteSF.h"
#include "MeshDataBaseInterface.h"

#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  enum evaluationType {
    UNKNOWNEVALTYPE,
    MEANRATIO,
    ORIENTEDMEANRATIO
  };

  // -------------------------------------------------------------------
  class elementEvaluatorBase {

  public:

    elementEvaluatorBase(const DiscreteSF * f): sizeField(f) {}
    virtual ~elementEvaluatorBase() {};
  
    void setSizeField(const DiscreteSF * sf) {sizeField = sf;}

    virtual int R_shape(const pRegion r, double * result) const = 0;
    virtual int R_shapeWithDisp(const pRegion rgn, double disp[4][3], 
                                double * shape) const = 0;
    virtual int F_shape(const pFace f, const double normal[3], double * result) const = 0;
    virtual int F_shapeWithDisp(const pFace pf, const double normal[3],
                                double disp[3][3], double * result) const = 0;
    virtual int XYZ_R_shape(const double[4][3], const pMSize[4], double * result) const = 0;
    virtual int XYZ_F_shape(const double[3][3], const pMSize[3], const double normal[3], double * result) const = 0;

    virtual inline double worstShapeEver() const {return 0.;}
    virtual inline double bestShapeEver()  const {return 1.;}
    virtual inline double whatsWorst(double shp1, double shp2) const {return std::min(shp1,shp2);}

    virtual evaluationType getType () const = 0;
    virtual std::string getName () const = 0;

  protected:
 
    const DiscreteSF * sizeField;
  };

  // -------------------------------------------------------------------

}

#endif
