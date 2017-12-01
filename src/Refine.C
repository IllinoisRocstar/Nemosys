#include <Refine.H>

Refine::Refine(meshBase* _mesh, std::string method, 
               int arrayID, double dev_mult, bool maxIsmin)
{
  mesh = _mesh;
  if (!mesh->getSFBool())
  {
    mesh->generateSizeField(method, arrayID, dev_mult, maxIsmin);
  }

  if (mesh->getSFBool())
  {
    mesh->writeMSH("converted.msh");
    vtkCellData* cd = mesh->getDataSet()->GetCellData();
    int i;
    if (cd)
    {
      std::string array_name = mesh->getDataSet()->GetPointData()->GetArrayName(arrayID);
      if (!method.compare("gradient"))
        array_name.append("GradientSF");
      else if (!method.compare("value"))
        array_name.append("ValueSF");

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
  else
  {
    std::cout << "Unable to generate size field for refinement" << std::endl;
    exit(1);
  }
}

Refine::~Refine()
{

  if (adapter)
  { 
    delete adapter;
    adapter = 0;
  }
  if (MadMesh)
  {
    MAd::M_delete(MadMesh);
    MadMesh = 0;
  }
  if (bSF) 
  { 
    delete bSF; 
    bSF = 0; 
  }
  remove("converted.msh");
  remove("backgroundSF.msh");
  remove("refined.msh");
  std::cout << "Refined destroyed" << std::endl;
}

void Refine::run()
{
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
  mesh->transfer(refinedVTK,"FE");
  std::string newname = mesh->getFileName();
  refinedVTK->write(trim_fname(newname, "_refined.vtu"), ".vtu");
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

