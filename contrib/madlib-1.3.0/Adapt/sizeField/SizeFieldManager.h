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

#ifndef _H_SIZEFIELDMANAGER
#define _H_SIZEFIELDMANAGER

#include "SizeFieldBase.h"
#include "PWLinearSField.h"
#include "LocalSizeField.h"
#include "AnalyticalSField.h"
#include "BackgroundSF.h"
#include "DistanceFunction.h"

#include <set>

namespace MAd {

  // -------------------------------------------------------------------
  // This class lists, intersects and updates all size fields
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------

  class SizeFieldManager {

  public:

    SizeFieldManager(pMesh m, pSField psf=NULL);
    ~SizeFieldManager();

  public:

    void setMesh(pMesh m);
    void setSmoothing(bool enable, double maxGrad) 
    {
      smooth = enable;
      maxGradient = maxGrad;
    }

    void intersect(const pSField sf);

    void addSizeField(pSField sf);

    void update();

    DiscreteSF* getSizeField()             { return mainSF; }
    const DiscreteSF* getSizeField() const { return mainSF; }

    std::set<LocalSizeField *> getLocalSizeFields() const { return locals; }

  private:

    void addPWLinearSF(PWLSField * sf);
    void addAnalyticalSF(AnalyticalSField * sf);
    void addLocalSF(LocalSizeField * sf);
    void addBGSF(BackgroundSF * bg);

  private:

    DiscreteSF* mainSF;

    std::set<PWLSField *>        linears;
    std::set<AnalyticalSField *> analyticals;
    std::set<LocalSizeField *>   locals;
    std::set<BackgroundSF *>       bgs;

    pMesh mesh;

    double localTime;

    bool smooth;
    double maxGradient;
  };

  // -------------------------------------------------------------------

}

#endif
