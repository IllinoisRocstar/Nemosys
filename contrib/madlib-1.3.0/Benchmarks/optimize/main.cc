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
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <stdlib.h>

using namespace MAd;

// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // Check input
  // ------------
  if ( argc != 2 ) {
    printf("Error: usage: \'executable mshFile\'\n");
    exit(0);
  }
  string meshFile = argv[1];

  // Build tools
  // ------------
  pGModel model = 0;
  pMesh mesh = M_new(model);
  M_load(mesh,meshFile.c_str());

  PWLSField * sizeField = new PWLSField(mesh);
  sizeField->setCurrentSize();
//   vector<string> h, e0, e1, e2;
//   h.push_back("0.50"); h.push_back("0.05"); h.push_back("0.01");
//   e0.push_back("1.0"); e0.push_back("0.0"); e0.push_back("0.0");
//   e1.push_back("0.0"); e1.push_back("1.0"); e1.push_back("0.0");
//   e2.push_back("0.0"); e2.push_back("0.0"); e2.push_back("1.0");
//   AnalyticalSField* asf = new AnalyticalSField(h, e0, e1, e2);
  
//   MeshAdapter* ma = new MeshAdapter(mesh,asf);
  MeshAdapter* ma = new MeshAdapter(mesh,sizeField);

  // Output situation before optimization
  // -------------------------------------
  printf ("Statistics before optimization: \n");
  ma->printStatistics(std::cout);
  ma->writePos("meanRatioBefore.pos",OD_MEANRATIO);

  // Optimize
  // ---------
  printf ("Optimizing mesh %s...\n\n",meshFile.c_str());
  ma->run();

  // Outputs final mesh
  // -------------------
  printf ("Statistics after optimization: \n");
  ma->printStatistics(std::cout);
  ma->writePos("meanRatioAfter.pos",OD_MEANRATIO);
  M_writeMsh (mesh, "result.msh", 2, NULL);

  // Cleaning
  // ---------
  if (ma) delete ma;
  if (sizeField) delete sizeField;
}

// ----------------------------------------------------------------------
