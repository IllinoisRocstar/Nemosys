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

#ifndef _H_MESSAGE_MAD_
#define _H_MESSAGE_MAD_

#include "MAdSingleton.h"

#include <ostream>
#include <list>

namespace MAd {

  // -------------------------------------------------------------------
  // abort function
  typedef void (*AbortFunction)(void *);

  // -------------------------------------------------------------------
  class MAdMsg {

  public:

    MAdMsg();
    ~MAdMsg() {}
  
    void initialize();
    void finalize();

    void registerAbortFct(AbortFunction, void *);

    void info   (int, const char*, const char*,...) const;
    void warning(int, const char*, const char*,...) const;
    void error  (int, const char*, const char*,...) const;

  private:

    std::string writePosition(int, const char*) const;

  private:

    std::ostream* outStream;
    std::ostream* errStream;

    std::list<std::pair<AbortFunction,void*> > abortFcts; 

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<MAdMsg> MAdMsgSgl;

}

#endif
