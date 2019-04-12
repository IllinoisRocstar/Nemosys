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

#ifndef _H_CALLBACKMANAGER
#define _H_CALLBACKMANAGER

#include "MAdSingleton.h"
#include "CallbackDefinition.h"

#include <list>

namespace MAd {

  // -------------------------------------------------------------------
  class CallBackManager {

  public:
  
    CallBackManager() {};
    ~CallBackManager() {};
  
    void initialize();
    void finalize();

    void registerCallBack(CBFunction CB, void* userData);
    void registerCallBackMove(CBFunc_move CB_move, void* userData_move);

    void unregisterCallBack(CBFunction CB, void* userData);
    void unregisterCallBackMove(CBFunc_move CB_move, void* userData_move);

    void callCallBacks(pPList before, pPList after,
                       operationType type , pEntity ppp);
    void callCallBackMoves(pVertex, double *);
  
  private:

    std::list<CBStructure> CB;
    std::list<CBStructureMove> CBMove;

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<CallBackManager> CallBackManagerSgl;

}

#endif
