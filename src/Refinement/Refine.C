#include "Refinement/Refine.H"

#include <vtkCellData.h>
#include <vtkPointData.h>

#include "AuxiliaryFunctions.H"
#include "Drivers/TransferDriver.H"
#include "Mesh/gmshMesh.H"

namespace NEM {
namespace ADP {

Refine::Refine(meshBase *_mesh, const std::string &method, int arrayID,
               double dev_mult, bool maxIsmin, double edge_scale,
               const std::string &_ofname, double sizeFactor, int order)
    : mesh(_mesh), ofname(_ofname), bndrConst(false) {
  if (method != "uniform") {
    // creates sizeField
    mesh->generateSizeField(method, arrayID, dev_mult, maxIsmin, sizeFactor,
                            order);
    initAdaptive(arrayID, method);
  } else if (method == "uniform") {
    initUniform(edge_scale);
  }
}

Refine::~Refine() {
  // destructor for adapter also destroys size field
  if (adapter) {
    delete adapter;
    adapter = nullptr;
  }
  /*
  if (bSF)
  {
    delete bSF;
    bSF = nullptr;
    std::cout << "here1" << std::endl;
  }
  if (pwlSF)
  {
    delete pwlSF;
    pwlSF = nullptr;
    std::cout << "here2" << std::endl;
  }
  */
  if (MadMesh) {
    MAd::M_delete(MadMesh);
    MadMesh = nullptr;
  }

  std::cout << "Refine destroyed" << std::endl;
}

void Refine::initUniform(double edge_scale) {
  std::cout << "Uniform Refinement Selected" << std::endl;

  auto *gMesh = new gmshMesh(mesh);
  gMesh->write("converted.msh", 2.0, false);

  MAd::pGModel gmodel = nullptr;

  MadMesh = MAd::M_new(gmodel);
  MAd::M_load(MadMesh, "converted.msh");
  classifyBoundaries();
  // DISCRETE/PWLSF SIZEFIELD
  pwlSF = new MAd::PWLSField(MadMesh);
  // sets size as mean edge length squared for edges adjacent to given vertex
  pwlSF->setCurrentSize();
  pwlSF->scale(edge_scale);
  // timing adapter construction
  nemAux::Timer T;
  T.start();
  // instantiating adapter with linear sf
  adapter = new MAd::MeshAdapter(MadMesh, pwlSF);

  // TJW
  // option to constrain boundary edges when refining
  if (bndrConst) adapter->setConstraint(bnd);

  T.stop();
  std::cout << "Time for adapter construction (ms): " << T.elapsed() << "\n \n";
  std::cout << "Refine constructed" << std::endl;
}

void Refine::initAdaptive(int arrayID, const std::string &method) {
  pwlSF = nullptr;
  auto *gMesh = new gmshMesh(mesh);
  gMesh->write("converted.msh", 2.0, false);

  vtkCellData *cd = mesh->getDataSet()->GetCellData();
  int i = 0;
  std::string array_name;
  if (cd) {
    array_name = mesh->getDataSet()->GetPointData()->GetArrayName(arrayID);
    if (method == "gradient")
      array_name.append("GradientSF");
    else if (method == "value")
      array_name.append("ValueSF");
    else if (method == "Z2 Error Estimator")
      array_name.append("Z2ErrorSF");
    for (i = 0; i < cd->GetNumberOfArrays(); ++i) {
      if (array_name == cd->GetArrayName(i)) break;
    }
    if (i == cd->GetNumberOfArrays()) {
      std::cerr << "Error: Did not find " << array_name << " in cell data set"
                << std::endl;
      exit(1);
    }
  }
  mesh->writeMSH("backgroundSF.msh", "cell", i, true);
  mesh->unsetCellDataArray(array_name);

  MAd::pGModel gmodel = nullptr;
  MadMesh = MAd::M_new(gmodel);
  MAd::M_load(MadMesh, "converted.msh");
  classifyBoundaries();

  bSF = new MAd::BackgroundSF("backgroundSF");
  bSF->loadData("backgroundSF.msh");

  std::cout << "\n \n Beginning Adapter Construction" << std::endl;
  // timing adapter construction
  nemAux::Timer T;
  T.start();
  // instantiating adapter with background sizeField
  adapter = new MAd::MeshAdapter(MadMesh, bSF);
  T.stop();
  std::cout << "Time for adapter construction (ms): " << T.elapsed() << "\n \n";
  std::cout << "Refine constructed" << std::endl;
}

void Refine::run(bool transferData, bool bndryConstraint) {
  if (!adapter) {
    std::cerr << "Adapter hasn't been constructed!" << std::endl;
    exit(1);
  }
  // set boundary constraint flag
  bndrConst = bndryConstraint;
  // option to constrain boundary edges when refining
  if (bndrConst) adapter->setConstraint(bnd);
  // Output situation before refinement
  std::cout << "Statistics before refinement: " << std::endl;
  adapter->printStatistics(std::cout);

  // Adaptive refinement
  std::cout << "Refining the mesh ..." << std::endl;
  adapter->run();
  std::cout << "Statistics after refinement: " << std::endl;
  adapter->printStatistics(std::cout);

  // running Laplacian smoothing
  std::cout << "Optimizing the mesh" << std::endl;

  adapter->checkTheMesh(10);

  for (int i = 0; i < 1; ++i) {
    adapter->LaplaceSmoothing();
    adapter->splitLongestEdges();
    adapter->removeSlivers();
    adapter->optimiseEdgeLength();
    adapter->optimiseElementShape();
  }
  /* there is a fast laplace smoothing function available where instead of
   * computing the optimal position it uses the cavity center. The center is the
   * initial position passed to the routine "computeOptimalLocation" before it
   * calculates optimal. The optimal way is run by default. If we modify
   * line 714/718 in AdaptInterface.cc to laplOp->runFast we can use the
   * fast method */

  // Outputs final mesh
  std::cout << "Statistics after optimization: " << std::endl;
  adapter->printStatistics(std::cout);

  // unclassifying boundary elements for proper output
  unClassifyBoundaries();
  // writing refined mesh to file in msh format
  MAd::M_writeMsh(MadMesh, "refined.msh", 2);
  meshBase *refinedVTK = meshBase::exportGmshToVtk("refined.msh");

  // mesh->setCheckQuality(1);

  if (transferData) {
    // mesh->transfer(refinedVTK, "Consistent Interpolation");
    auto transfer = NEM::DRV::TransferDriver::CreateTransferObject(
        mesh, refinedVTK, "Consistent Interpolation");
    transfer->run(mesh->getNewArrayNames());
  }

  refinedVTK->setFileName(ofname);
  refinedVTK->report();
  refinedVTK->write();

  delete refinedVTK;
}

void Refine::classifyBoundaries() {
  // finding and classifying boundary elements as 2d
  MadMesh->classify_unclassified_entities();
  MadMesh->destroyStandAloneEntities();
  bnd = (MAd::pGEntity)MAd::GM_faceByTag(MadMesh->model, 0);
  MadMesh->classify_grid_boundaries(bnd);
}

void Refine::unClassifyBoundaries() {
  // unclassifying boundary elements for proper output
  MadMesh->unclassify_grid_boundaries();
}

}  // namespace ADP
}  // namespace NEM
