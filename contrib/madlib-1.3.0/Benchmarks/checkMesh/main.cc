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
#include <sys/resource.h>
#include <string>
using std::string;
#include <stdlib.h>

using namespace MAd;

#ifdef PARALLEL
#include "mpi.h"
#endif

// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {

#ifdef PARALLEL 
  MPI_Init(&argc, &argv);
#endif

  // Check input
  // ------------
  if ( argc != 3 && argc != 2  ) {
    printf("Error: usage: \'executable mshFile [geoFile]\'\n");
    exit(0);
  }
  string meshFile = argv[1];

  // Build tools
  // ------------
  printf ("Checking mesh %s...\n\n",meshFile.c_str());

  // --- Reading model ---
  pGModel model = 0;
  GM_create(&model,"theModel");
  
  if ( argc == 3 ) {
    string geoFile  = argv[2];
    GM_read(model, geoFile.c_str());
  }
  else {
    GM_readFromMSH(model, meshFile.c_str());
  }

  pMesh mesh = M_new(model);
  M_load(mesh,meshFile.c_str());

  MeshStatus status = VALID;
  if ( !checkMesh(mesh, CHECK_ALL, 1, std::cout, &status) )
    {
      printf("\nThe mesh is invalid, status = %d\n\n",status);
    }
  else
    {
      printf("\nThe mesh is valid, status = %d\n\n",status);
    }


#ifdef PARALLEL 
  MPI_Finalize();
#endif
}

// ----------------------------------------------------------------------
