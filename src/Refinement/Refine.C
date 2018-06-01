#include <Refine.H>
#include <vtkPointData.h>
#include <vtkCellData.h>

Refine::Refine(meshBase* _mesh, const std::string& method, 
               int arrayID, double dev_mult, bool maxIsmin, 
               double edge_scale, std::string _ofname, double sizeFactor)
{
  mesh = _mesh;
  ofname = _ofname; 
  if (!mesh->getSFBool() && method.compare("uniform"))
  {
    // creates sizefield and switches sfbool
    mesh->generateSizeField(method, arrayID, dev_mult, maxIsmin, sizeFactor);
  }
  if (!method.compare("uniform"))
  {
    initUniform(edge_scale);
  } 
  if (mesh->getSFBool())
  {
    initAdaptive(arrayID,method);
  }
  else if (method.compare("uniform"))
  {
    std::cout << "Unable to generate size field for refinement" << std::endl;
    exit(1);
  }
}

Refine::~Refine()
{
  // destructor for adapter also destroys size field
  if (adapter)
  { 
    delete adapter;
    adapter = 0;
  }
/*  if (bSF) 
  { 
    delete bSF; 
    bSF = 0;
    std::cout << "here1" << std::endl; 
  }
  if (pwlSF)
  {
    delete pwlSF;
    pwlSF = 0;
    std::cout << "here2" << std::endl; 
  }*/

  if (MadMesh)
  {
    MAd::M_delete(MadMesh);
    MadMesh = 0;
  }
  remove("converted.msh");
  remove("backgroundSF.msh");
  remove("refined.msh");
  std::cout << "Refine destroyed" << std::endl;
}

void Refine::initUniform(double edge_scale)
{
  std::cout << "Uniform Refinement Selected" << std::endl;
  mesh->writeMSH("converted.msh");
  MAd::pGModel gmodel = 0;
  MadMesh = MAd::M_new(gmodel);
  MAd::M_load(MadMesh,"converted.msh");
  classifyBoundaries();
  // DISCRETE/PWLSF SIZEFIELD
  pwlSF = new MAd::PWLSField(MadMesh);
  // sets size as mean edge length squared for edges adjacent
  // to given vertex  
  pwlSF->setCurrentSize();
  pwlSF->scale(edge_scale);
  // timing adapter construction
  Timer T;
  T.start();    
  // instantiating adapter with linear sf
  adapter = new MAd::MeshAdapter(MadMesh, pwlSF);
  T.stop();
  std::cout << "Time for adapter construction (ms): " << T.elapsed() << "\n \n";
  std::cout << "Refine constructed" << std::endl;
  // for safety
  mesh->setSFBool(0);
}

void Refine::initAdaptive(int arrayID, const std::string& method)
{
  pwlSF = 0;
  mesh->writeMSH("converted.msh");
  vtkCellData* cd = mesh->getDataSet()->GetCellData();
  int i;
  std::string array_name;
  if (cd)
  {
    array_name = mesh->getDataSet()->GetPointData()->GetArrayName(arrayID);
    if (!method.compare("gradient"))
      array_name.append("GradientSF");
    else if (!method.compare("value"))
      array_name.append("ValueSF");
    else if (!method.compare("Z2 Error Estimator"))
      array_name.append("Z2ErrorSF");
    for (i = 0; i < cd->GetNumberOfArrays(); ++i)
    {
      std::string currname = cd->GetArrayName(i);
      if (!array_name.compare(currname))
        break;
    }
    if (i == cd->GetNumberOfArrays())
    {
      std::cout << "Error: Did not find " << array_name << " in cell data set" << std::endl;
      exit(1);
    }
  }
  mesh->writeMSH("backgroundSF.msh", "cell", i, 1);
  mesh->unsetCellDataArray(&array_name[0u]); 

  MAd::pGModel gmodel = 0;
  MadMesh = MAd::M_new(gmodel);
  MAd::M_load(MadMesh,"converted.msh");
  classifyBoundaries();
       
  bSF = new MAd::BackgroundSF("backgroundSF");
  bSF->loadData("backgroundSF.msh"); 
  
  std::cout << "\n \n Beginning Adapter Construction" << std::endl;
  // timing adapter construction
  Timer T;
  T.start();
  // instantiating adapter with background sizefield
  adapter = new MAd::MeshAdapter(MadMesh, bSF);
  T.stop();
  std::cout << "Time for adapter construction (ms): " << T.elapsed() << "\n \n";
  std::cout << "Refine constructed" << std::endl;
}

void Refine::run(bool transferData)
{
  if (!adapter)
  {
    std::cout << "Adapter hasn't been constructed!" << std::endl;
    exit(1);
  }
  // Output situation before refinement
  std::cout << "Statistics before refinement: " << std::endl;
  adapter->printStatistics(std::cout);

  // Adaptive refinement
  std::cout << "Refining the mesh ..." << std::endl;
  adapter->run();
  std::cout << "Statistics after refinement: " << std::endl;
  adapter->printStatistics(std::cout);

  // running laplacian smoothing
  std::cout << "Optimizing the mesh" << std::endl;

  for (int i = 0; i < 1; ++i)
  {
    adapter->LaplaceSmoothing();
    adapter->splitLongestEdges();
    adapter->removeSlivers();
    adapter->optimiseEdgeLength();
    adapter->optimiseElementShape();    
  }
   /* there is a fast laplace smoothing function available where instead of computing 
     the optimal position it uses the cavity center. The center is the initial position
     passed to the routine "computeOptimalLocation" before it calculates optimal. The 
     optimal way is run by default. If we modify line 714/718 in AdaptInterface.cc to
     laplOp->runFast we can use the fast method */  


  // Outputs final mesh
  std::cout << "Statistics after optimization: " << std::endl;
  adapter->printStatistics(std::cout);
  
  // unclassifying boundary elements for proper output
  unClassifyBoundaries();
  // writing refined mesh to file in msh format 
  MAd::M_writeMsh(MadMesh, "refined.msh", 2);
  
  meshBase* refinedVTK = meshBase::exportGmshToVtk("refined.msh");
  //mesh->setCheckQuality(1);
  if (transferData)
    mesh->transfer(refinedVTK,"Consistent Interpolation");

  refinedVTK->setFileName(ofname);
  refinedVTK->report();
  refinedVTK->write();
   
  if (refinedVTK)
  {
    delete refinedVTK;
    refinedVTK = 0; 
  }
}

void Refine::classifyBoundaries()
{
  // finding and classifying boundary elements as 2d 
  MadMesh->classify_unclassified_entities();
  MadMesh->destroyStandAloneEntities();
  MAd::pGEntity bnd = (MAd::pGEntity) MAd::GM_faceByTag(MadMesh->model, 0);
  MadMesh->classify_grid_boundaries(bnd);
}

void Refine::unClassifyBoundaries()
{
  // unclassifying boundary elements for proper output
  MadMesh->unclassify_grid_boundaries();
}

