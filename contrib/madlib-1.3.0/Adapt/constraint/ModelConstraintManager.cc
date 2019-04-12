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

#include "ModelConstraintManager.h"

namespace MAd {

  // -------------------------------------------------------------------
  void ModelConstraintManager::initialize(pGModel _model)
  {
    model = _model;
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::finalize()
  {
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::setModel(pGModel _model)
  {
    model = _model;
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::clear()
  {
    constrEntities.clear();
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::constrain(int type, int id)
  {
    pGEntity pGE = GM_entityByTag(model, type, id);
    constrain(pGE);
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::constrain(pGEntity pGE)
  {
    // constrain entity
    constrEntities.insert(pGE);
  
#ifdef _HAVE_GMSH_
    // constrain lower geometrical levels
    std::list<pGEntity> subGE = GEN_closure(pGE);
    std::list<pGEntity>::const_iterator subIter = subGE.begin();
    for (; subIter != subGE.end(); subIter++) constrain(*subIter);
#endif
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::unconstrain(int type, int id)
  {
    pGEntity pGE = GM_entityByTag(model, type, id);
    unconstrain(pGE);
  }

  // -------------------------------------------------------------------
  void ModelConstraintManager::unconstrain(pGEntity pGE)
  {
    std::set<pGEntity>::iterator eIter = constrEntities.find(pGE);
    if ( eIter != constrEntities.end() ) constrEntities.erase(eIter);
  }

  // -------------------------------------------------------------------
  bool ModelConstraintManager::constrained(pGEntity pGE)
  {
    //GCTODO: check upper-level geometric entities: need model with connectivity

    if ( constrEntities.find(pGE) != constrEntities.end() ) return true;
    return false;
  }

  // -------------------------------------------------------------------

}
