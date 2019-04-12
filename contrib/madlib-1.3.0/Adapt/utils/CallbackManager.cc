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

#include "CallbackManager.h"

using std::list;

namespace MAd {

  // -------------------------------------------------------------------
  void CallBackManager::initialize()
  {
  }

  // -------------------------------------------------------------------
  void CallBackManager::finalize()
  {
    CB.clear();
    CBMove.clear();
  }

  // -------------------------------------------------------------------
  void CallBackManager::registerCallBack(CBFunction CBFct, 
                                         void* userData)
  {
    CBStructure cb;
    cb.function = CBFct;
    cb.data = userData;
    CB.push_back(cb);
  }

  // -------------------------------------------------------------------
  void CallBackManager::registerCallBackMove(CBFunc_move CB_move, 
                                             void* userData_move)
  {
    CBStructureMove cbm;
    cbm.function = CB_move;
    cbm.data = userData_move;
    CBMove.push_back(cbm);
  }

  // -------------------------------------------------------------------
  void CallBackManager::unregisterCallBack(CBFunction CBFct, 
                                           void* userData)
  {
    CBStructure cb;
    cb.function = CBFct;
    cb.data = userData;
    CB.remove(cb);
  }

  // -------------------------------------------------------------------
  void CallBackManager::unregisterCallBackMove(CBFunc_move CB_move, 
                                               void* userData_move)
  {
    CBStructureMove cbm;
    cbm.function = CB_move;
    cbm.data = userData_move;
    CBMove.remove(cbm);
  }

  // -------------------------------------------------------------------
  void CallBackManager::callCallBacks(pPList before, pPList after,
                                      operationType type , pEntity ppp)
  {
    list<CBStructure>::const_iterator cbIter = CB.begin();
    list<CBStructure>::const_iterator cbLast = CB.end();

    for (; cbIter != cbLast; cbIter++) {
      CBFunction f = (*cbIter).function;
      void * userData = (*cbIter).data;
      f(before,after,userData,type,ppp);
    }
  }

  // -------------------------------------------------------------------
  void CallBackManager::callCallBackMoves(pVertex pv, double * xyz)
  {
    list<CBStructureMove>::const_iterator cbIterM = CBMove.begin();
    list<CBStructureMove>::const_iterator cbLastM = CBMove.end();

    for (; cbIterM != cbLastM; cbIterM++) {
      CBFunc_move fM = (*cbIterM).function;
      void * userDataM = (*cbIterM).data;
      fM(pv,xyz,userDataM);
    }
  }

  // -------------------------------------------------------------------

}
