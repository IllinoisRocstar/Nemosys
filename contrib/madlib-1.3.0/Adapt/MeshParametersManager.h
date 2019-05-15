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

#ifndef _H_MESHPARAMETERSMANAGER
#define _H_MESHPARAMETERSMANAGER

#include "MeshDataBaseInterface.h"
#include "MAdOperatorBase.h"
#include "MAdSingleton.h"

namespace MAd {

  // -------------------------------------------------------------------
  class MeshParametersManager {
    
  public:
  
    MeshParametersManager();
    ~MeshParametersManager() {};

    void initialize();
    void finalize();

    // ------------------------------------------
  public:
  
    void setLowerLengthSqBound(double); 
    void setUpperLengthSqBound(double); 
    void setSliverTriBound(double);
    void setSliverTetBound(double);
  
    double getLowerLengthSqBound() const {return lowerLengthSqBound;}
    double getUpperLengthSqBound() const {return upperLengthSqBound;}
    double getSliverBound(int dim) const;
    double getSliverBound(const pMesh) const;
    double getSliverTriBound() const {return sliverTriQualityBound;}
    double getSliverTetBound() const {return sliverTetQualityBound;}

  private:

    double lowerLengthSqBound, upperLengthSqBound;
    double sliverTetQualityBound;
    double sliverTriQualityBound;

    // ------------------------------------------
  public:

    void setSliverPermissionInESplit    (bool, double bound=-1.);
    void setSliverPermissionInECollapse (bool, double bound=-1.);

    bool getSliverPermissionInESplit    () const {return sliverPermInESplit;}
    bool getSliverPermissionInECollapse () const {return sliverPermInECollapse;}

    double getSliverLowerLengthSqBound  () const {return sliverLowerLengthSqBound;}
    double getSliverUpperLengthSqBound  () const {return sliverUpperLengthSqBound;}

  private:

    bool sliverPermInESplit, sliverPermInECollapse;
    double sliverLowerLengthSqBound, sliverUpperLengthSqBound;

    // ------------------------------------------
  public:

    void setNoSwapQuality(double noSwap) {noSwapQuality = noSwap;}
    void setSwapMinImproveRatio(double ratio) {swapMinImproveRatio = ratio;}

    double getNoSwapQuality() const {return noSwapQuality;}
    double getSwapMinImproveRatio() const {return swapMinImproveRatio;}

  private:

    double noSwapQuality;
    double swapMinImproveRatio;

    // ------------------------------------------
  public:
    
    void setBigLength(double len) { bigLength = len; }
    double getBigLength() { return bigLength; }

  private:
    
    double bigLength;

    // ------------------------------------------
  public:

    void diagnostics() const;

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<MeshParametersManager> MeshParametersManagerSgl;

}

#endif
