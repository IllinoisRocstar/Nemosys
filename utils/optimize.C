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

// Gmsh headers 
#include <GModel.h>
#include "CGNSOptions.h"

// MAdLib headers 
#include "MAdLib.h"
#include "NodalDataManager.h"

// Nemosys headers
#include <vtkAnalyzer.H> 

// others
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace MAd;

// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {

  GModel* srcGModel;
  srcGModel = new GModel("source"); 

  // Check input
  if ( argc != 2 ) {
    printf("Error: usage: \'executable vtkFile\'\n");
    exit(0);
  }
  std::string srcVtkFile = argv[1];

  // read input vtk file and convert to gmsh mesh file
  srcGModel->readVTK(srcVtkFile);
  srcGModel->writeMSH("vtkConverted.msh", 2.2);
  
  // write statistics of input mesh
  //srcGModel->writePOS("meshData.pos", true, true, true, true, true, true, true);

  // read physical quantities in mesh
  vtkAnalyzer* srcVTK;
  srcVTK = new vtkAnalyzer(argv[1]);
  srcVTK->read();
  srcVTK->report();

  // Build tools for mesh adaptation (AMR)
  pGModel model = 0;
  pMesh mesh = M_new(model);
  M_load(mesh,"vtkConverted.msh");

  // inqure mesh object
  std::cout << "---------- Input Mesh Information --------------" << std::endl;
  std::cout << "Number of points in the mesh    = " << mesh->nbPoints << std::endl;
  std::cout << "Number of Edges in the mesh     = " << mesh->nbEdges << std::endl;
  std::cout << "Number of Triangles in the mesh = " << mesh->nbTriangles << std::endl;
  std::cout << "Number of Quads in the mesh     = " << mesh->nbQuads << std::endl;
  std::cout << "Number of Tets in the mesh      = " << mesh->nbTets << std::endl;
  std::cout << "Number of Hexes in the mesh     = " << mesh->nbHexes << std::endl;
  std::cout << "Number of Prisms in the mesh    = " << mesh->nbPrisms << std::endl;
  std::cout << "------------------------------------------------" << std::endl;

  /*
  PWLSField * sizeField = new PWLSField(mesh);
  //sizeField->setCurrentSize();
  sizeField->setAllVSizes(0.1);
  MeshAdapter* ma = new MeshAdapter(mesh,sizeField);
  */

  
  std::vector<std::string> h, e0, e1, e2;
  h.push_back("0.05+0.1*y"); h.push_back("0.1"); h.push_back("0.1");
  e0.push_back("1.0"); e0.push_back("0.0"); e0.push_back("0.0");
  e1.push_back("0.0"); e1.push_back("1.0"); e1.push_back("0.0");
  e2.push_back("0.0"); e2.push_back("0.0"); e2.push_back("1.0");
  //AnalyticalSField* sizeField = new AnalyticalSField(h, e0, e1, e2);
  AnalyticalSField* sizeField = 
          new AnalyticalSField("-0.18*sin(x*3.14/2)*sin(y*3.14)+0.2");
  MeshAdapter* ma = new MeshAdapter(mesh,sizeField);
  

  // attach some data to mesh
  std::vector<double> pntData;
  pntData.resize(srcVTK->getNumberOfPoints(), 0.0);
  int pntId = 0;
  for(std::vector<double>::iterator it = pntData.begin(); it != pntData.end(); ++it){
   double* pntCoords = srcVTK->getPointCoords(pntId++); 
    //*it = pntCoords[0]*pntCoords[0];
    *it = pntCoords[0];
  }
  ma->registerData("T", pntData);
  NodalDataManagerSgl::instance().writeData("T","solutionPre.pos");

  // make sure data are attached properly
  std::vector<double> preAMRData;
  ma->getMeshData("T", &preAMRData);
  std::cout << "preAMRData[0] = " << preAMRData[0] << std::endl;
  std::cout << "Size of vector = " << preAMRData.size() << std::endl; 
  double summation = 0.0;
  for(std::vector<double>::iterator it = preAMRData.begin(); it != preAMRData.end(); ++it)
    summation += *it;
  std::cout << "sumVector = " << summation << std::endl;
  std::cout << "Average value = " << summation*1.0 / preAMRData.size() << std::endl;

  // write attached dataset to a vtk file
  srcVTK->setPointDataArray("T", 1, preAMRData);
  srcVTK->report();
  srcVTK->write("oldSolution.vtu");
  
  // write old solution to csv file
  //srcVTK->writeCSV("oldSolution.csv", preAMRData);
  
  // using MAd to write mesh and data into CSV file
  NodalDataManagerSgl::instance().writeDataCSV("T", "oldSolutionMAd.csv");

  // Output situation before optimization
  printf ("Statistics before optimization: \n");
  ma->printStatistics(std::cout);
  //ma->writePos("meanRatioBefore.pos",OD_MEANRATIO);

  // adjust adaptor settings
  ma->setSwapOnBoundary();
  //ma->setSliverQuality(0.5);

  // Optimize
  printf ("Optimizing the mesh ...\n");
  ma->run();
  ma->removeSlivers();

  // Outputs final mesh
  printf ("Statistics after optimization: \n");
  ma->printStatistics(std::cout);
  //ma->writePos("meanRatioAfter.pos",OD_MEANRATIO);
  M_writeMsh (mesh, "newMesh.msh", 2, NULL);
  
  // inqure mesh object
  std::cout << "---------- Output Mesh Information -------------" << std::endl;
  std::cout << "Number of points in the mesh    = " << mesh->nbPoints << std::endl;
  std::cout << "Number of Edges in the mesh     = " << mesh->nbEdges << std::endl;
  std::cout << "Number of Triangles in the mesh = " << mesh->nbTriangles << std::endl;
  std::cout << "Number of Quads in the mesh     = " << mesh->nbQuads << std::endl;
  std::cout << "Number of Tets in the mesh      = " << mesh->nbTets << std::endl;
  std::cout << "Number of Hexes in the mesh     = " << mesh->nbHexes << std::endl;
  std::cout << "Number of Prisms in the mesh    = " << mesh->nbPrisms << std::endl;
  std::cout << "------------------------------------------------" << std::endl;

  // get data after refinement
  std::vector<double> postAMRData;
  ma->getMeshData("T", &postAMRData);
  std::cout << "postAMRData[0] = " << postAMRData[0] << std::endl;
  std::cout << "Size of vector = " << postAMRData.size() << std::endl; 
  summation = 0.0;
  for(std::vector<double>::iterator it = postAMRData.begin(); it != postAMRData.end(); ++it)
    summation += *it;
  std::cout << "sumVector = " << summation << std::endl;
  std::cout << "Average value = " << summation*1.0/postAMRData.size() << std::endl;

  // write mesh with field data to GMSH POS files
  NodalDataManagerSgl::instance().writeData("T","solutionPost.pos");

  // using MAd to write mesh and data into CSV file
  NodalDataManagerSgl::instance().writeDataCSV("T", "newSolutionMAd.csv");

  // convert final Gmsh file back to vtk
  GModel* trgGModel;
  trgGModel = new GModel("targer"); 
  trgGModel->readMSH("newMesh.msh");
  trgGModel->writeVTK("newMesh.vtk", false, true);

  // write physical quantities to vtk file
  vtkAnalyzer* trgVTK;
  trgVTK = new vtkAnalyzer((char*)"newMesh.vtk");
  trgVTK->read();
  trgVTK->report();
  trgVTK->setPointDataArray("T", 1, postAMRData);
  trgVTK->report();
  trgVTK->write("newSolution.vtu");

  // Cleaning
  if (ma) delete ma;
  if (sizeField) delete sizeField;
  if (srcVTK) delete srcVTK;
  if (srcGModel) delete srcGModel;
  if (trgGModel) delete trgGModel;
  if (trgVTK) delete trgVTK;

}

// ----------------------------------------------------------------------
