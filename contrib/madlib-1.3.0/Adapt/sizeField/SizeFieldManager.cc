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

#include "SizeFieldManager.h"
#include "NullSField.h"
#include "MAdTimeManager.h"
#include "MeshQualityManager.h"
#include "MAdMessage.h"
#include "MeshSizeBase.h"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::set;

namespace MAd {

  // -------------------------------------------------------------------
  SizeFieldManager::SizeFieldManager(pMesh m, pSField psf): 
    localTime(-1.),smooth(false),maxGradient(1.)
  {
    if (!m) return;
    mesh = m;

    mainSF = new PWLSField(mesh,"Main size field");

    localTime = MAdTimeManagerSgl::instance().getTime();

    if (!psf) MAdMsgSgl::instance().warning(__LINE__,__FILE__,"No size field registered");

    else {
    
      mainSF->intersect(psf);
  
      sFieldType type = psf->getType();
      switch (type) {
      case NULLSFIELD: {
        break;
      }
      case DISCRETESFIELD: {
        assert ( ((DiscreteSF*)psf)->discretization() == VERTEX_P1_DSFTYPE );
        PWLSField* linsf = static_cast<PWLSField*>(psf);
        addPWLinearSF( linsf );
        break;
      }
      case ANALYTICALSFIELD: {
        AnalyticalSField* asf = static_cast<AnalyticalSField*>(psf);
        addAnalyticalSF( asf );
        break;
      }
      case LOCALSFIELD: {
        LocalSizeField* locsf = static_cast<LocalSizeField*>(psf);
        addLocalSF( locsf );
        break;
      }
      case BACKGROUNDSFIELD: {
        BackgroundSF* bg = static_cast<BackgroundSF*>(psf);
        addBGSF( bg );
        break;
      }
      default: {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "unknown size field type %d", type);
      }
      }
    }  
  }

  // -------------------------------------------------------------------
  SizeFieldManager::~SizeFieldManager() {

    if (mainSF) delete mainSF;
    linears.clear();
    analyticals.clear();
    locals.clear();
    bgs.clear();
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::setMesh(pMesh m)
  {
    mesh = m;
    mainSF = new PWLSField(mesh);
    MAdMsgSgl::instance().warning(__LINE__,__FILE__,"No size field registered");
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::intersect(const pSField sf) {
    mainSF->intersect(sf);
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::addSizeField(pSField sf) {

    mainSF->intersect(sf);
  
    sFieldType type = sf->getType();
    switch (type) {
    case NULLSFIELD: {
      break;
    }
    case DISCRETESFIELD: {
      assert ( ((DiscreteSF*)sf)->discretization() == VERTEX_P1_DSFTYPE );
      PWLSField* linsf = static_cast<PWLSField*>(sf);
      addPWLinearSF( linsf );
      break;
    }
    case ANALYTICALSFIELD: {
      AnalyticalSField* asf = static_cast<AnalyticalSField*>(sf);
      addAnalyticalSF( asf );
      break;
    }
    case LOCALSFIELD: {
      LocalSizeField* locsf = static_cast<LocalSizeField*>(sf);
      addLocalSF( locsf );
      break;
    }
    case BACKGROUNDSFIELD: {
      BackgroundSF* bg = static_cast<BackgroundSF*>(sf);
      addBGSF( bg );
      break;
    }
    default: {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Unknown size field type %d",type);
    }
    }
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::addPWLinearSF(PWLSField* sf) {
    linears.insert(sf);
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::addAnalyticalSF(AnalyticalSField* sf) {
    analyticals.insert(sf);
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::addLocalSF(LocalSizeField* sf) {
    locals.insert(sf);
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::addBGSF(BackgroundSF* sf) {
    bgs.insert(sf);
  }

  // -------------------------------------------------------------------
  void SizeFieldManager::update()
  {
    MeshQualityManagerSgl::instance().clearAllShapes(); // not useful in isotropic case (be careful about the hypothesis)

    mainSF->cleanUp();

    set<PWLSField*>::const_iterator linIter = linears.begin();
    set<PWLSField*>::const_iterator linLast = linears.end();
    for (; linIter != linLast; linIter++) {
      mainSF->intersect(*linIter);
    }

    set<AnalyticalSField*>::const_iterator aIter = analyticals.begin();
    set<AnalyticalSField*>::const_iterator aLast = analyticals.end();
    for (; aIter != aLast; aIter++) {
      mainSF->intersect(*aIter);
    }

    set<LocalSizeField*>::iterator locIter = locals.begin();
    set<LocalSizeField*>::iterator locLast = locals.end();
    for (; locIter != locLast; locIter++) {
      (*locIter)->updateTree();
      mainSF->intersect(*locIter);
    }

    set<BackgroundSF*>::iterator bgIter = bgs.begin();
    set<BackgroundSF*>::iterator bgLast = bgs.end();
    for (; bgIter != bgLast; bgIter++) {
      mainSF->intersect(*bgIter);
    }

    if ( smooth ) mainSF->smooth(maxGradient);

    localTime = MAdTimeManagerSgl::instance().getTime();
  }

  // -------------------------------------------------------------------

}
