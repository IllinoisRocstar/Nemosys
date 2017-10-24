#include <MAdLib.h>
#include <GModel.h>
#include <vtkAnalyzer.H>
#include <meshPhys.H>
#include <NodalDataManager.h>
#include <chrono>

int main(int argc, char* argv[])
{
  // Check input
  if ( argc != 5 ) {
    std::cout << "Usage Error: " << argv[0] << " meshFile outputMesh array_id size" << std::endl;
    exit(1);
  }
  std::string meshFile = argv[1];


  // ------------------- MESHPHYS TESTS ----------------------------//
  std::cout << "TESTING MESHPHYS CLASS" << std::endl; 
  meshPhys* mshphys;
  mshphys = new meshPhys((char*) &(meshFile)[0u]);
  PointDataArray& pntdata = mshphys->getPointData(0);
 
  /*std::cout << "Name of arrays " << pntdata.getName() << std::endl; 
  std::cout << "NumComponents: " <<  pntdata.getNumComponent() << std::endl;
  std::cout << "NumTuples: " << pntdata.getNumTuple() << std::endl;
  std::cout << "first few values:" << std::endl;
  for (int i = 0; i < pntdata.getNumTuple(); ++i)
  { 
    for (int j = 0; j < pntdata.getNumComponent(); ++j)
    {
      std::cout << mshphys->getPointData(0)(i,j) << " ";
    }
    std::cout << std::endl;
    if (i > 50)
      break;
  } */ 

  int array_id = std::stoi(argv[3]);
  double size = std::stof(argv[4]); 
 
  std::vector<double> gradients = mshphys->ComputeL2GradAtAllCells(array_id);
  reciprocal_vec(gradients); 
  std::vector<double> lengths = mshphys->GetCellLengths();
  std::vector<double> gradminmax = getMinMax(gradients);
  std::vector<double> lengthminmax = getMinMax(lengths);
  std::cout << "Gradient and Length Tests: \n";
  for (int i = 0; i < mshphys->numberOfCells; ++i)
  {   
    std::cout << lengths[i] << " " << gradients[i] <<  std::endl;
  }

  lengthminmax[1] -= (lengthminmax[1] - lengthminmax[0])/2.0;
  std::cout << "\n \n Minmax Lengths " << lengthminmax[0] << " "
            << lengthminmax[1] << std::endl;  

  std::cout << "\n \n Minmax Gradients " << gradminmax[0] << " "
            << gradminmax[1] << std::endl;  

  scale_vec_to_range(gradients, gradminmax, lengthminmax); 
  std::vector<double> scalesminmax = getMinMax(gradients);
  std::cout << "Scale Min and Max: ";
  std::cout << scalesminmax[0] << " " << scalesminmax[1] << std::endl;
  for (int i = 0; i < mshphys->numberOfCells; ++i)
    std::cout << gradients[i] << std::endl;  
  

  mshphys->getCellsWithPoint(10);
  
  mshphys->report(); 
 // delete mshphys;
  //--------------------------------- END MESHPHYS TESTS ----------------------//

  // create gmodel from input mesh
  MAd::pGModel gmodel = 0;
  MAd::pMesh mesh = MAd::M_new(gmodel);
  vtkAnalyzer* tmpmesh;
  if (meshFile.find(".v") != -1 ) 
  {
    tmpmesh = new vtkAnalyzer((char*) &(meshFile)[0u]);
    tmpmesh->read();
    std::string ofname = "converted.msh";
    tmpmesh->writeMSH(ofname);

    MAd::M_load(mesh, "converted.msh"); 
  }
  else if (meshFile.find(".msh") != -1)
  {
    MAd::M_load(mesh, (char*) &(meshFile)[0u]); 
  }

  // find 2d boundary elements
  MAd::pGEntity bnd = (MAd::pGEntity) MAd::GM_faceByTag(mesh->model, 0);
  mesh->classify_grid_boundaries(bnd);
  mesh->classify_unclassified_entities();
  mesh->destroyStandAloneEntities();

  // writing background mesh and taking care of non tet/tri entities
  mshphys->writeBackgroundMSH("backgroundSF.msh", gradients);
  delete mshphys;
  MAd::BackgroundSF* bSF = new MAd::BackgroundSF("backgroundSF");
  bSF->loadData("backgroundSF.msh");
  std::cout << "\n Beginning Adapter Construction \n" << std::endl;
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, bSF);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for adapter construction (ms): " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() 
            << std::endl;
  
  /* TODO: writing vector data not implemented in madlib
           instead, split vector data into component vectors
           if writing data required */  


  // retrieving point data from vtk mesh
  std::vector<std::vector<double>> pntData;
  int numTuple, numComponent; 
  tmpmesh->getPointDataArray(array_id, pntData, numTuple, numComponent);

  // attaching point data to adapter's mesh
  adapter->registerVData("displacement", pntData);
  // check if data is correctly attached
  std::vector<std::vector<double>> preAMRData;
  adapter->getMeshVData("displacement", &preAMRData);
  std::cout << "preAMRData[1][1] = " << preAMRData[1][1] << std::endl;
  std::cout << "Size of vector = " << preAMRData.size() << std::endl; 

  // Output situation before optimization
  std::cout << "Statistics before optimization: " << std::endl;
  adapter->printStatistics(std::cout);
  
  // Optimize
  std::cout << "Optimizing the mesh ..." << std::endl;
  adapter->run();
  adapter->removeSlivers();

  // Outputs final mesh
  std::cout << "Statistics after optimization: " << std::endl;
  adapter->printStatistics(std::cout);

  // writing refined mesh to file in msh format 
  MAd::M_writeMsh(mesh, argv[2], 2);
   
  // get data after refinement
  std::vector<std::vector<double>> postAMRData;
  adapter->getMeshVData("displacement", &postAMRData);
  std::cout << "postAMRData[1][1] = " << postAMRData[1][1] << std::endl;
  std::cout << "Size of vector = " << postAMRData.size() << std::endl; 

  // convert refined msh file back to vtk
  GModel* trgGModel;
  trgGModel = new GModel("refined"); 
  trgGModel->readMSH("refined.msh");
  trgGModel->writeVTK("refined.vtk", false, true); // binary=false, saveall=true

  // write physical quantities to vtk file
  vtkAnalyzer* trgVTK;
  trgVTK = new vtkAnalyzer((char*) "refined.vtk");
  trgVTK->read();
  trgVTK->report();
  std::vector<double> postAMRDatas = flatten(postAMRData);

  trgVTK->setPointDataArray("displacement", 3, postAMRDatas);
  trgVTK->report();
  trgVTK->write("refined_solution.vtu");
  trgVTK->writeMSH("refined_solution.msh");

 
  //if (localSF) delete localSF;
  //if (adapter) delete adapter;
  if (bSF) delete bSF;
  if (trgGModel) delete trgGModel;
  if (trgVTK) delete trgVTK;
  return 0;
}


  /*MAd::PWLSField* pwlSF = new MAd::PWLSField(mesh);
  pwlSF->setCurrentSize();
  pwlSF->scale(.5);
  MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, pwlSF);
  */
  

  /*size_t beg = 0;
  size_t end = meshFile.find(".");
  string name; 
  if (end != -1) 
  { 
    name = meshFile.substr(beg,end);
    name.append(".vtu");
    std::cout << "writing vtu ... " << std::endl;
    tmpmesh->write((char*) &(name)[0u]);
    std::cout << "done" << std::endl;
    
  }
  else 
  {
    std::cout << "Error finding file extension for " << meshFile << std::endl;
    exit(1);
  }*/

/*   // Reading gmsh file and info
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

  srcGModel->readVTK(meshFile);
  std::cout << "Num Vertices: " << srcGModel->getNumVertices() << std::endl;

  std::set<GVertex*, GEntityLessThan>::iterator iter = srcGModel->firstVertex();  
  */

  /*/////////////////// LOCAL SIZEFIELD TESTS /////////////////////////    
  MAd::LocalSizeField* localSF = new MAd::LocalSizeField(mesh);
  //localSF->setIsoSize(0.1,".01");
  localSF->setIsoSize(0.01,"0.01");
  localSF->addGeometricEntity(3,1);
  localSF->addGeometricEntity(3,2);
  localSF->addGeometricEntity(3,3);
  localSF->addGeometricEntity(3,4);
  localSF->addGeometricEntity(3,5);
  localSF->updateTree();
  */

  // FOR 2D
  /*localSF->setIsoSize(0.0001,"0.0");
  localSF->addGeometricEntity(2,33);
  localSF->addGeometricEntity(2,34);
  localSF->addGeometricEntity(2,35);
  localSF->addGeometricEntity(2,36);
  localSF->addGeometricEntity(2,37);
  localSF->updateTree();
  */

  /*MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, localSF);


  // DISCRETE/PWLSF SIZEFIELD
  // Must define if local doesn't cover full domain
  MAd::PWLSField* pwlSF = new MAd::PWLSField(mesh);
  // sets size as mean edge length squared for edges adjacent
  // to given vertex  
  pwlSF->setCurrentSize();
  //pwlSF->setAllVSizes(1.0);
  pwlSF->scale(1);
  
  //MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, pwlSF);

  adapter->addSizeField(pwlSF);
  */
  // FOR 2D
  /*MAd::LocalSizeField* localSF2 = new MAd::LocalSizeField(mesh);
  localSF2->setIsoSize(0.1,"0.01");
  localSF2->addGeometricEntity(2,33);
  localSF2->addGeometricEntity(2,34);
  localSF2->addGeometricEntity(2,35);
  localSF2->addGeometricEntity(2,36);
  localSF2->addGeometricEntity(2,37);
  localSF2->updateTree();
 
  adapter->addSizeField(localSF2);*/ // was used to refine around inclusion, but leave inclusion empty 
  ////////////////////////// END LOCAL SIZEFIELD TESTS ////////////////////
