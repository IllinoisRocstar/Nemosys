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

#ifdef _HAVE_GMSH_
#include "gmsh/GmshGlobal.h"
#endif

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#endif

#include "MAdMessage.h"
#include "MAdResourceManager.h"

#include <string>
using std::string;

using namespace MAd;

void MAdLibInitialize(int * argc, char** argv[])
{
#ifdef PARALLEL
  MPI_Init(argc,argv);
  AP_init(argc,argv);
#endif

#ifdef _HAVE_GMSH_
  char * tmp[1];
  tmp[0] = (char*)" ";
  GmshInitialize(1,tmp);
  unsigned int i = 1;
  int j = 0;
  GmshSetOption("General","Terminal",i,j);
#endif

  MAdMsgSgl::instance().initialize();
  MAdResourceManagerSgl::instance().initialize();
}

void MAdLibFinalize()
{
  MAdResourceManagerSgl::instance().finalize();
  MAdMsgSgl::instance().finalize();

#ifdef _HAVE_GMSH_
  GmshFinalize();
#endif

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  AP_finalize();
  MPI_Finalize();
#endif
}

