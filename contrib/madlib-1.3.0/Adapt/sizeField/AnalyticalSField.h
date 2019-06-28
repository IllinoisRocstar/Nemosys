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

#ifndef _H_ANALYTICALSFIELD
#define _H_ANALYTICALSFIELD

#include "madlib_export.h"

#include "SizeFieldBase.h"

#include <vector>

namespace MAd {

  class MAdStringFieldEvaluator;

  // -------------------------------------------------------------------
  typedef pMSize (*sizeFunction)(const double[3],double);

  // -------------------------------------------------------------------
  class MADLIB_EXPORT AnalyticalSField: public SizeFieldBase
  {
  public:

    AnalyticalSField();
    AnalyticalSField(std::string);
    AnalyticalSField(std::vector<std::string>,std::vector<std::string>,
                     std::vector<std::string>,std::vector<std::string>);
    AnalyticalSField(sizeFunction);
    ~AnalyticalSField();

  public:

    sFieldType getType() const { return ANALYTICALSFIELD; }
    void describe() const;

    void scale(double);

    // set the size
    void setSize(sizeFunction);
    void setSize(const std::string);
    void setSize(std::vector<std::string>,std::vector<std::string>,
                 std::vector<std::string>,std::vector<std::string>);

    // get the size at a location (allocate space!)
    pMSize getSize(const pVertex) const;
    pMSize getSizeOnEntity(const pEntity, const double[3]) const;

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

  private:

    // description of the sizes if we use functions
    sizeFunction sFct;

    // description of the sizes if we use strings
    bool isotropic;
    std::string h0, h1, h2;
    std::vector<std::string> e0, e1, e2;
    MAdStringFieldEvaluator * evaluator;

  private:

    pMSize eval(const double[3]) const;
  };

}

// -------------------------------------------------------------------
#endif
