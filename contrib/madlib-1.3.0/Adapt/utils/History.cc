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

#include "History.h"

#include <stdarg.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using std::ostream;

namespace MAd {

  static int actionIntParameters[UNKNOWN_ACTION] = {1,1,1,1,1};

  // -------------------------------------------------------------------
  void History::initialize()
  {
    open = false;
    refCount = -1;
  }

  // -------------------------------------------------------------------
  void History::finalize()
  {
    journal.clear();
    refJournal.clear();
  }

  // -------------------------------------------------------------------
  void History::add(int ot, actionType at, ...)
  {
    if ( !open ) return;

    va_list ap;
    va_start(ap, at);

    std::vector<int> intParams;
    for (int iInt=0; iInt < actionIntParameters[at]; iInt++) {
      intParams.push_back(va_arg(ap, int));
    }

    va_end(ap);

    action act;
    act.opType    = ot;
    act.actType   = at;
    act.intParams = intParams;

    add(act);
  }

  // -------------------------------------------------------------------
  void History::add(action act)
  {
    if ( !open ) return;

    journal.push_back(act);
  
    if ( !compare() ) {
      std::cerr << "Error in History: operations are not the same as in the reference journal\n";
      flushJournal("journal_for_diff");
      exit(1);
    }
  }

  // -------------------------------------------------------------------
  bool History::compare()
  {
    if ( !open ) return false;

    if ( refJournal.empty() || refCount >= (int) refJournal.size() ) return true;
  
    for (; refCount < (int) journal.size(); refCount++) {
      if ( journal[refCount] != refJournal[refCount] ) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------
  void History::loadJournal(std::string name)
  {
    FILE * f = fopen(name.c_str(),"r");
    if ( !f ) {
      std::cerr << "Error: could not open file " << name << std::endl; throw;
    }

    int opType;
    int actType;
    while ( fscanf(f,"%d %d",&opType,&actType) != EOF ) {

      std::vector<int> intParams;
      int par;
      for (int iPar=0; iPar<actionIntParameters[actType]; iPar++) {
        fscanf(f,"%d",&par);
        intParams.push_back(par);
      }

      action act;
      act.opType    = opType;
      act.actType   = (actionType)    actType;
      act.intParams = intParams;

      refJournal.push_back(act);
    }

    fclose(f);

    refCount = 0;
  }

  // -------------------------------------------------------------------
  void History::flushJournal(ostream& out) const
  {
    for (int iAct=0; iAct < (int) journal.size(); iAct++) {
      action act = journal[iAct];
      out << act.opType << "\t" << act.actType;
      for (int iP=0; iP<actionIntParameters[act.actType]; iP++) {
        out << "\t"<< act.intParams[iP];
      }
      out << "\n";
    }
  }

  // -------------------------------------------------------------------
  void History::flushJournal(std::string name) const
  {
    std::ofstream fileOut(name.c_str());
    flushJournal(fileOut);
  }

  // -------------------------------------------------------------------

}
