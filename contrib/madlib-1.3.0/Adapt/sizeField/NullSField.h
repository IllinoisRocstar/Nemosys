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

#ifndef _H_NULLSFIELD
#define _H_NULLSFIELD

#include "SizeFieldBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class NullSField : public SizeFieldBase
  {
  public:

    NullSField(): SizeFieldBase() {}
    ~NullSField() {}

    sFieldType getType() const { return NULLSFIELD; }

    // edge length (squared)
    double SF_VV_lengthSq(const pVertex, const pVertex) const;
    double SF_XYZ_lengthSq(const double[3], const double[3],
                           const pMSize, const pMSize=NULL) const;

    // face area (squared)
    double SF_F_areaSq(const pFace) const;
    double SF_XYZ_areaSq(const double[3][3], const pMSize,
                         const double[3]) const;

    // region volume
    double SF_R_volume(const pRegion) const;
    double SF_XYZ_volume(const double[4][3], const pMSize) const;

    // center and its associated size
    double SF_E_center(const pEdge, double[3], double * reducSq, pMSize *) const;
    double SF_VV_center(const pVertex, const pVertex,
                        double[3], double * reducSq, pMSize *) const;

  public: // functions that have no sense here...

    void scale(double) {return;}

    void setSize(pEntity, pMSize) {return;}
    void setSize(pEntity, double) {return;}
    void deleteSize(pEntity) {return;}

    pMSize getSize(const pVertex) const {return NULL;}
    pMSize getSizeOnEntity(const pEntity, const double[3]) const {return NULL;};
  };

}

// -------------------------------------------------------------------

#endif
