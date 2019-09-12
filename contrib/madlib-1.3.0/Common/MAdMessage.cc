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

#include "MAdMessage.h"

#include <cstdarg>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>

using std::stringstream;
using std::string;
using std::ostream;

namespace MAd {

  // -------------------------------------------------------------------
  MAdMsg::MAdMsg()
  {
    outStream = &std::cout;
    errStream = &std::cerr;
  }

  // -------------------------------------------------------------------
  void MAdMsg::initialize()
  {
  }

  // -------------------------------------------------------------------
  void MAdMsg::finalize()
  {
  }

  // -------------------------------------------------------------------
  void MAdMsg::registerAbortFct(AbortFunction fct, void * data)
  {
    abortFcts.push_back(std::make_pair(fct,data));
  }

  // -------------------------------------------------------------------
  void MAdMsg::info(int line,
                    const char* file,
                    const char* fmt,...) const
  {
    char buff[1024];
    va_list args;
    va_start (args, fmt);
    vsprintf(buff, fmt, args);
    va_end (args);

    if ( line >= 0 ) {
      *outStream << "Info: " << buff << writePosition(line,file) << "" << std::endl;
    }
    else {
      *outStream << "Info: " << buff << std::endl;
    }
  }

  // -------------------------------------------------------------------
  void MAdMsg::warning(int line,
                       const char* file,
                       const char* fmt,...) const
  {
    char buff[1024];
    va_list args;
    va_start (args, fmt);
    vsprintf(buff, fmt, args);
    va_end (args);

    if ( line >= 0 ) {
      *outStream << "WARNING: " << buff << writePosition(line,file) << "" << std::endl;
    }
    else {
      *outStream << "WARNING: " << buff << std::endl;
    }
  }

  // -------------------------------------------------------------------
  void MAdMsg::error(int line,
                     const char* file,
                     const char* fmt,...) const
  {
    char buff[1024];
    va_list args;
    va_start (args, fmt);
    vsprintf(buff, fmt, args);
    va_end (args);

    if ( line >= 0 ) {
      *outStream << "ERROR: " << buff << writePosition(line,file) << "" << std::endl;
    }
    else {
      *outStream << "ERROR: " << buff << std::endl;
    }

    // call abort functions
    std::list<std::pair<AbortFunction,void*> >::const_iterator iter = abortFcts.begin();
    for (; iter != abortFcts.end(); iter++) {
      AbortFunction fct = (*iter).first;
      void * data       = (*iter).second;
      fct(data);
    }

    abort();
  }

  // -------------------------------------------------------------------
  string MAdMsg::writePosition(int line, const char* file) const
  {
    stringstream ss;
    string iStr;  ss << line;  ss >> iStr;
    return  " (line " + iStr + " in file \'" + file + "\')";
  }

  // -------------------------------------------------------------------

}
