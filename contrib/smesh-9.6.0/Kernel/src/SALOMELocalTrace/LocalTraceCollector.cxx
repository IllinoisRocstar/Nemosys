// Copyright (C) 2007-2020  CEA/DEN, EDF R&D, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

//  File   : LocalTraceCollector.cxx
//  Author : Paul RASCLE (EDF)
//  Module : KERNEL
//  $Header$
//
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include "LocalTraceCollector.hxx"

// ============================================================================
/*!
 *  This class is for use without CORBA, inside or outside SALOME.
 *  SALOME uses SALOMETraceCollector, to allow trace collection via CORBA.
 *  Type of trace (and corresponding class) is chosen in LocalTraceBufferPool.
 *
 *  Guarantees a unique object instance of the class (singleton thread safe)
 *  a separate thread for loop to print traces is launched.
 */
// ============================================================================

BaseTraceCollector* LocalTraceCollector::instance()
{
  if (_singleton == 0) // no need of lock when singleton already exists
    {
      _singletonMutex.lock(); // acquire lock to be alone
      if (_singleton == 0)                     // another thread may have got
        {                                      // the lock after the first test
          BaseTraceCollector* myInstance = new LocalTraceCollector();

          _thread = new std::thread(LocalTraceCollector::run, nullptr);
          _sem.wait();
          _singleton = myInstance; // _singleton known only when init done
        }
      _singletonMutex.unlock(); // release lock
    }
  return _singleton;
}

// ============================================================================
/*!
 *  In a separate thread, loop to print traces.
 *  Mutex guarantees initialisation on instance method is done and only one run
 *  allowed (double check ...)
 *  Loop until there is no more buffer to print,
 *  and no ask for end from destructor.
 *  Get a buffer. If type = ABORT then exit application with message.
 */
// ============================================================================

void* LocalTraceCollector::run(void *bid)
{
  _sem.post(); // unlock instance

  LocalTraceBufferPool* myTraceBuffer = LocalTraceBufferPool::instance();
  LocalTrace_TraceInfo myTrace;

  // --- Loop until there is no more buffer to print,
  //     and no ask for end from destructor.

  while ((!_threadToClose) || myTraceBuffer->toCollect() )
    {
      if (_threadToClose)
        {
          DEVTRACE("FileTraceCollector _threadToClose");
          //break;
        }

      myTraceBuffer->retrieve(myTrace);
      if (myTrace.traceType == ABORT_MESS)
        {
          std::cout << std::flush ;
#ifndef WIN32
          std::cerr << "INTERRUPTION from thread " << myTrace.threadId
               << " : " <<  myTrace.trace;
#else
          std::cerr << "INTERRUPTION from thread " << (void*)(&myTrace.threadId)
               << " : " <<  myTrace.trace;
#endif
          std::cerr << std::flush ; 
          exit(1);     
        }
      else if (myTrace.traceType == NORMAL_MESS)
        {
          std::cout << std::flush ;
#ifndef WIN32
          std::cerr << "th. " << myTrace.threadId << " " << myTrace.trace;
#else
          std::cerr << "th. " << (void*)(&myTrace.threadId)
               << " " << myTrace.trace;
#endif
          std::cerr << std::flush ; 
        }
      else 
        {
          std::cout << std::flush ;
#ifndef WIN32
          std::cerr << myTrace.trace;
#else
          std::cerr << myTrace.trace;
#endif
          std::cerr << std::flush ; 
        }
    }
  return NULL;
}

// ============================================================================
/*!
 *  Destructor: wait until printing thread ends (LocalTraceCollector::run)
 */
// ============================================================================

LocalTraceCollector:: ~LocalTraceCollector()
{
  _singletonMutex.lock(); // acquire lock to be alone
  if (_singleton)
    {
      DEVTRACE("LocalTraceCollector:: ~LocalTraceCollector()");
      LocalTraceBufferPool* myTraceBuffer = LocalTraceBufferPool::instance();
      _threadToClose = 1;
      myTraceBuffer->insert(NORMAL_MESS,"end of trace\n"); // to wake up thread
      if (_thread)
        {
          _thread->join();
          delete _thread;
          _thread = 0;
          _threadToClose = 0;
        }
      _singleton = 0;
    }
  _singletonMutex.unlock(); // release lock
}

// ============================================================================
/*!
 * Constructor: no need of LocalTraceBufferPool object initialization here,
 * thread safe singleton used in LocalTraceBufferPool::instance()
 */
// ============================================================================

LocalTraceCollector::LocalTraceCollector()
{
  _thread = 0;
  _threadToClose = 0;
}


