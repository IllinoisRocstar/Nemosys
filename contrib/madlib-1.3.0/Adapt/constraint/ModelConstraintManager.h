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

#ifndef _H_MODELCONSTRAINTMANAGER
#define _H_MODELCONSTRAINTMANAGER

#include "ModelInterface.h"
#include <set>

#include "MAdSingleton.h"

namespace MAd {

  // -------------------------------------------------------------------

  class ModelConstraintManager {

  public:
  
    ModelConstraintManager()  {};
    ~ModelConstraintManager() {};

    void initialize(pGModel);
    void finalize();
    void setModel(pGModel);

  public:

    void clear();

    void constrain  (int type, int id);
    void constrain  (pGEntity);

    void unconstrain(int type, int id);
    void unconstrain(pGEntity);

    bool constrained(pGEntity);

  private:

    pGModel model;
    std::set<pGEntity> constrEntities;
  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<ModelConstraintManager> ModelConstraintManagerSgl;

}

#endif
