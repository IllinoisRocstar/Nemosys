#include <MAdLib.h>
#include <GModel.h>
#include <vtkAnalyzer.H>

int main(int argc, char* argv[])
{
  /*MAd::pGModel gmodel = NULL;
  MAd::GM_create(&gmodel, "GeoModel");
  MAd::GM_read(gmodel, argv[1]);
  
  MAd::pMesh mesh = MAd::M_new(gmodel);
  MAd::M_load(mesh, argv[2]); 
  
  MAd::PWLSField* sizeField = new MAd::PWLSField(mesh);
  
  sizeField->setCurrentSize();
  sizeField->scale(0.5);
  
  MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, sizeField);
  adapter->run();

  MAd::M_writeMsh(mesh, "newMesh.msh", 2);
  */

   


  // Check input
  if ( argc != 2 ) {
    std::cout << "Usage Error: " << argv[0] << "vtkFile" << std::endl;
    exit(1);
  }
  std::string VTKFile = argv[1];
/*
   // Reading gmsh file and info
    std::string mshFileName = argv[6];
    std::ifstream mshStream(mshFileName.c_str());

    MAd::pGModel model = NULL;
    MAd::GM_create(&model,"initial");
    MAd::pMesh mesh = MAd::M_new(model);
    MAd::M_load(mesh, mshFileName.c_str());
    std::cout << "Using tags?: " << MAd::GM_physical(model) << std::endl;
  
  //  std::cout << model.getNumVertices() << std::endl;

    std::vector<std::vector<double>> curds;
    MAd::VIter vit = M_vertexIter(mesh);
    MAd::pVertex pv;
    while (pv = VIter_next(vit)) {
      MAd::pPoint pp = MAd::V_point(pv);
      std::cout << pp->X << " " << pp->Y << " " << pp->Z << std::endl;
    }*/

  // read input vtk file and convert to gmsh mesh file
  /*GModel* srcGModel;
  srcGModel = new GModel("source"); 

  srcGModel->readVTK(VTKFile);
  std::cout << "Num Vertices: " << srcGModel->getNumVertices() << std::endl;

  std::set<GVertex*, GEntityLessThan>::iterator iter = srcGModel->firstVertex();  
  */

  vtkAnalyzer* tmpmesh;
  tmpmesh = new vtkAnalyzer((char*) &(VTKFile)[0u]);
  tmpmesh->read();
 
  /*size_t beg = 0;
  size_t end = VTKFile.find(".");
  string name; 
  if (end != -1) 
  { 
    name = VTKFile.substr(beg,end);
    name.append(".vtu");
    std::cout << "writing vtu ... " << std::endl;
    tmpmesh->write((char*) &(name)[0u]);
    std::cout << "done" << std::endl;
    
  }
  else 
  {
    std::cout << "Error finding file extension for " << VTKFile << std::endl;
    exit(1);
  }

  srcGModel->writeMSH("vtkConverted.msh", 2.2);*/

  std::string ofname = "testmsh.msh";
  tmpmesh->writeMSH(ofname);
 
  return 0;
}



