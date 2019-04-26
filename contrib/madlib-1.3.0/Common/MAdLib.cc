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
  #include <gmsh.h>
#endif

#ifdef PARALLEL
  #include "mpi.h"
  #include "autopack.h"
#endif

#include "MAdMessage.h"
#include "MAdResourceManager.h"

#include <string>

void MAdLibInitialize(int *argc, char **argv[]) {
#ifdef PARALLEL
  MPI_Init(argc,argv);
  AP_init(argc,argv);
#endif

#ifdef _HAVE_GMSH_
  char *tmp[1];
  tmp[0] = (char *) " ";
  gmsh::initialize(1, tmp, false);
  unsigned int i = 1;
  gmsh::option::setNumber("General.Terminal", i);
#endif

  MAd::MAdMsgSgl::instance().initialize();
  MAd::MAdResourceManagerSgl::instance().initialize();
}

void MAdLibFinalize() {
  MAd::MAdResourceManagerSgl::instance().finalize();
  MAd::MAdMsgSgl::instance().finalize();

#ifdef _HAVE_GMSH_
  gmsh::finalize();
#endif

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  AP_finalize();
  MPI_Finalize();
#endif
}
