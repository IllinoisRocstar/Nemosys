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

#ifndef _H_HISTORY
#define _H_HISTORY

#include "MAdSingleton.h"

#include <vector>
#include <ostream>
#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  enum actionType {
    OP_CHECKCONSTRAINTS = 0,
    OP_CHECKGEOMETRY    = 1,
    OP_CHECKSHAPES      = 2,
    OPERATOR_APPLY      = 3,
    OPERATOR_UNAPPLY    = 4,
    UNKNOWN_ACTION      = 5
  };

  struct action {
    int opType;
    actionType actType;
    std::vector<int> intParams;

    bool operator!=(const action& other) {
      if ( opType  != other.opType  ) return true;
      if ( actType != other.actType ) return true;
      if ( intParams.size() != other.intParams.size() ) return true;
      for (unsigned int i=0; i<intParams.size(); i++) {
        if ( intParams[i] != other.intParams[i] ) return true;
      }
      return false;
    }
  };

  // -------------------------------------------------------------------
  class History {

  public:
  
    History()  {};
    ~History() {};
  
    void initialize();
    void finalize();
  
    void openJournal()  { open = true;  };
    void closeJournal() { open = false; };

    void add(int, actionType,...);
    void add(action);

    bool compare();

    void flushJournal(std::ostream& out) const;
    void flushJournal(std::string name) const;

    void loadJournal(std::string name);
  
  private:
  
    bool open; // indicates if the journal is in use

    std::vector<action> journal;      // list of actions

    std::vector<action> refJournal;  // journal to compare with
    int refCount;                    // position of the next action to compare
  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<History> HistorySgl;

}

#endif
