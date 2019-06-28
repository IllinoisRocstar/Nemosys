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

#include "MAdResourceManager.h"
#include "MAdMessage.h"

#include <fstream>

#ifdef PARALLEL
 #include <mpi.h>
#else
 #ifdef HAVE_UNISTD_H
  #include <sys/time.h>
  #include <unistd.h>
 #else
  #include <time.h>
  #include <windows.h> 
  #if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
   #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
  #else
   #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
  #endif

struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
int get_time_of_day(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
 
  if (NULL != tv)
    {
      GetSystemTimeAsFileTime(&ft);
 
      tmpres |= ft.dwHighDateTime;
      tmpres <<= 32;
      tmpres |= ft.dwLowDateTime;
 
      // converting file time
      tmpres -= DELTA_EPOCH_IN_MICROSECS; 
      tmpres /= 10;  // convert into microseconds
      tv->tv_sec = (long)(tmpres / 1000000UL);
      tv->tv_usec = (long)(tmpres % 1000000UL);
    }
 
  if (NULL != tz)
    {
      if (!tzflag)
        {
          _tzset();
          tzflag++;
        }
      tz->tz_minuteswest = _timezone / 60;
      tz->tz_dsttime = _daylight;
    }
 
  return 0;
}
 #endif

#endif

#if defined(__linux) || defined(linux)

void process_mem_usage(double& vm_usage, double& resident_set)
{
  using std::ios_base;
  using std::ifstream;
  using std::string;
  
  vm_usage     = 0.0;
  resident_set = 0.0;
  
  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);
  
  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;
  
  // the two fields we want
  //
  unsigned long vsize;
  long rss;
  
  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
  
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024.0 / 1024;
  resident_set = rss * page_size_kb / 1024;
}

#endif

// -------------------------------------------------------------------

namespace MAd {

  // -------------------------------------------------------------------
  void MAdResourceManager::initialize()
  {
  }

  // -------------------------------------------------------------------
  void MAdResourceManager::finalize()
  {
  }

  // -------------------------------------------------------------------
  double MAdResourceManager::getTime() const {
  
#ifdef PARALLEL
    return MPI_Wtime();
#else
    struct timeval tp;
    struct timezone tz;
#ifdef HAVE_UNISTD_H
    gettimeofday(&tp, &tz);
#else
    get_time_of_day(&tp, &tz);
#endif

    return ((double) tp.tv_sec +
            (double) ((double) .000001 * (double) tp.tv_usec));
#endif
  }

  // -------------------------------------------------------------------
  double MAdResourceManager::elapsedTime() const {

    return ( getTime() - initialTime );

  }

  // -------------------------------------------------------------------
  bool MAdResourceManager::getMemoryUsage(double& resident_size,
                                          double& virtual_size) const
  {
#if defined(__linux) || defined(linux)
    process_mem_usage(virtual_size,resident_size);
  
#ifdef PARALLEL
    double send[2] = {virtual_size,resident_size};
    double recv[2];
    MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    virtual_size  = recv[0];
    resident_size = recv[1];
#endif
    return true;
#else
    MAdMsgSgl::instance().info(__LINE__,__FILE__,
                               "Process memory statistics not reliable on this platform");
    return false;
#endif
  }

  // -------------------------------------------------------------------
  void MAdResourceManager::printMemoryUsage(std::string step, 
                                            std::ostream& out) const
  {
    double ram, virt;
    getMemoryUsage(ram,virt);
    out << "Memory usage at step \'"<<step<<"\': " 
        << ram  << " Mb (resident), "
        << virt << " Mb (virtual)\n";
  }
}

// -------------------------------------------------------------------
