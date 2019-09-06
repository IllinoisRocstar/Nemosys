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

#include "MAdLib.h"

#include <iostream>
using std::cout;
#include <sys/time.h>
#include <sys/resource.h>
#include <string>
using std::string;
#include <stdlib.h>
#include <fstream>

#if defined(__linux) || defined(linux)
  #include <unistd.h>
#endif

using namespace MAd;

#ifdef PARALLEL
#include "mpi.h"
#endif

// ----------------------------------------------------------------------
#ifdef PARALLEL
class EmptyExchanger : public MDB_DataExchanger
{
public:
  EmptyExchanger(int _tag): MDB_DataExchanger(_tag) {}
  ~EmptyExchanger() {}

  void * sendData (pEntity pe,    // in
		   int iProcDest, // in
		   int &_size ) {
    _size = 0;
    return NULL;
  }
  void receiveData (pEntity pe,      //in
	            int iProcSender, //in
		    void *buf ) {}
  void deleteExternalData( pEntity pe) const {}
};
#endif


// ----------------------------------------------------------------------
double CPUTime() {
  
#ifdef PARALLEL
  return MPI_Wtime();
#else
  struct timeval tp;
  struct timezone tz;

  gettimeofday(&tp,&tz);

  return ((double) tp.tv_sec +
          (double) ((double) .000001 * (double) tp.tv_usec));
#endif
}

// ----------------------------------------------------------------------
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

// ----------------------------------------------------------------------
void displayMemoryUsage(string step="")
{
  double ram, virt;

#if defined(__linux) || defined(linux)
  process_mem_usage(virt,ram);
  
#ifdef PARALLEL
    double send[2] = {virt,ram};
    double recv[2];
    MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    virt  = recv[0];
    ram = recv[1];
#endif

  std::cout << "Memory usage at step \'"<<step<<"\': " 
            << ram  << " Mb (resident), "
            << virt << " Mb (virtual)\n";
#endif
}

// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {

  MAdLibInitialize(&argc, &argv);

  displayMemoryUsage("initialized");

  // Check input
  // ------------
  if ( argc != 3 && argc != 2  ) {
    printf("Error: usage: \'executable mshFile [geoFile]\'\n");
    exit(0);
  }
  string meshFile = argv[1];

  // Build tools
  // ------------
  printf ("Analyzing mesh %s...\n\n",meshFile.c_str());

  // --- Reading model ---
  pGModel model = 0;
  GM_create(&model,"theModel");
  
  if ( argc == 3 ) {
    string geoFile  = argv[2];
    GM_readFromGEO(model, geoFile.c_str());
  }
  else {
    GM_readFromMSH(model, meshFile.c_str());
  }

  displayMemoryUsage("built model");

  pMesh mesh = M_new(model);
  M_load(mesh,meshFile.c_str());

  displayMemoryUsage("built mesh");

  SizeFieldBase * sizeField = new AnalyticalSField("1.");
  
  MeshAdapter* ma = new MeshAdapter(mesh,sizeField);

  displayMemoryUsage("built adapter");

  // temporary info
//   int nbVolEdges = 0;
//   int dim = M_dim(mesh);
//   int numR[20];
//   for (int i=0; i<20; i++) numR[i]=0;
//   EIter ei = M_edgeIter(mesh);
//   while ( pEdge edge = EIter_next(ei) ) {
//     if ( E_whatInType(edge) == dim ) {
//       int num = E_numRegions(edge);
//       numR[num] = numR[num]++;
//       nbVolEdges++;
//     }
//   }
//   EIter_delete(ei);
//   printf("Num points in edges crown:\n");
//   printf("Edges: %d\n",nbVolEdges);
  
//   for (int i=0; i<20; i++) {
//     double ratio = ( (double) (numR[i]) ) / ( (double)nbVolEdges );
//     printf("%d\t%d\t%f\n",i,numR[i],ratio);
//   }
//   printf("\n");


  // Outputs
  // --------
//   M_writeMsh(mesh,"test.msh",2);
  ma->printStatistics(std::cout);
  ma->writePos("meanRatio.pos",OD_MEANRATIO);

  // Cleaning
  // ---------
  if (ma) delete ma;
  if (sizeField) delete sizeField;
  M_delete(mesh);
  GM_delete(model);

  displayMemoryUsage("cleaned");

  MAdLibFinalize();
}

// ----------------------------------------------------------------------
