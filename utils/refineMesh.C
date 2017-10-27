#include <MAdLib.h>
#include <GModel.h>
#include <vtkAnalyzer.H>
#include <meshPhys.H>
#include <NodalDataManager.h>
#include <chrono>

// TODO: Investigate data transfer problem. Check out implementation in other utilities
//       Maybe you just fucked it up?

//---------------------Auxiliary Functions--------------------------------//
// computes reciprocal of argument
template <typename T>
T reciprocal(T x);

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x);

// find minmax of vector, excluding inf 
std::vector<double> getMinMax(const std::vector<double>& x);

// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
double scale_to_range(double x, const std::vector<double>& xminmax, 
                                const std::vector<double>& yminmax);

// scales each number in vector in range [xmin, xmax] to [ymin,ymax]
void scale_vec_to_range(std::vector<double>& x, 
                        const std::vector<double>& xminmax,
                        const std::vector<double>& yminmax);

// multiplies vector by scalar
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x);

// get average and stdev of values
std::vector<double> getMeanStdev(std::vector<double>& x);

// generates boolean array with 1 if value >= tol, 0 otherwise
std::vector<int> cellsToRefine(std::vector<double>& values, double tol);

// string trimming for consistent file names 
std::string trim_fname(std::string name, std::string ext);
//----------------------------------------------------------------------//


int main(int argc, char* argv[])
{
  // Check input
  if ( argc != 4 ) {
    std::cout << "Usage Error: " << argv[0] << " meshFile outputMesh array_id" << std::endl;
    exit(1);
  }
  std::string meshFile = argv[1];
  std::string outMeshFile = argv[2];
  
  // ------------------- MESHPHYS TESTS ----------------------------//
  std::cout << "TESTING MESHPHYS CLASS" << std::endl; 
  meshPhys* mshphys;
  mshphys = new meshPhys((char*) &(meshFile)[0u]);

  int array_id = std::stoi(argv[3]);
  std::string array_name =  mshphys->getPointData(array_id).getName(); 
  int dim = mshphys->getDimArray(array_id);
  
  std::vector<double> lengths = mshphys->GetCellLengths();
  std::vector<double> lengthminmax = getMinMax(lengths);

  lengthminmax[1] = lengthminmax[0];
  lengthminmax[0] = lengthminmax[0] - lengthminmax[0]/2;
  
  // for 1D data, use L2 to write abs
  std::vector<double> values = mshphys->ComputeL2ValAtAllCells(array_id);
  // get mean and stdev of values 
  std::vector<double> meanStdev = getMeanStdev(values);
  // get bool array of which cells to refine based on multiplier of stdev
  std::vector<int> cells2Refine = cellsToRefine(values, meanStdev[0]+meanStdev[1]/2.);
  // normalize values by mean
  std::vector<double> values_norm = (1./meanStdev[0])*values;  
  // take the reciprocal of values for size definition (high value -> smaller size)
  reciprocal_vec(values);
  // scale values to min max circumsphere diam of cells 
  // now, values represents a size field
  std::vector<double> valuesMinMax = getMinMax(values);
  scale_vec_to_range(values, valuesMinMax, lengthminmax);
  // setting sizes (values) to min element diam based on return of cellsToRefine function
  for (int i = 0; i < values.size(); ++i)
  {
    if (!cells2Refine[i])
      values[i] = lengthminmax[1];
  } 

 
  /*mshphys->setCellDataArray(&array_name[0u], dim, values); 
  mshphys->report();
  std::string newname;
  newname.insert(0,array_name);
  newname.insert(0, "_");
  newname.append("AtCell.vtu");
  mshphys->write((char*) &(trim_fname(meshFile, newname))[0u]);
 */
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
  mshphys->writeBackgroundMSH("backgroundSF.msh", values);
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

  if (dim > 1)
  {
    adapter->registerVData(array_name, pntData);
    // check if data is correctly attached
    std::vector<std::vector<double>> preAMRData;
    adapter->getMeshVData(array_name, &preAMRData);
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
    MAd::M_writeMsh(mesh, (char*) &outMeshFile[0u], 2);
   
    // get data after refinement
    std::vector<std::vector<double>> postAMRData;
    adapter->getMeshVData(array_name, &postAMRData);
    std::cout << "postAMRData[1][1] = " << postAMRData[1][1] << std::endl;
    std::cout << "Size of vector = " << postAMRData.size() << std::endl; 

    // convert refined msh file back to vtk
    GModel* trgGModel;
    trgGModel = new GModel("refined"); 
    trgGModel->readMSH((char*) &outMeshFile[0u]);
    trgGModel->writeVTK(trim_fname(outMeshFile,".vtk"), false, true); // binary=false, saveall=true


    // write physical quantities to vtk file
    vtkAnalyzer* trgVTK;
    trgVTK = new vtkAnalyzer((char*) &(trim_fname(outMeshFile,".vtk")[0u]));
    trgVTK->read();
    std::vector<double> postAMRDatas = flatten(postAMRData);
    trgVTK->setPointDataArray(&array_name[0u], dim, postAMRDatas);
    trgVTK->report();
    trgVTK->write((char*) &(trim_fname(outMeshFile, "_solution.vtu"))[0u]);
    trgVTK->writeMSH(trim_fname(outMeshFile, "_solution.msh"));
    if (trgGModel) delete trgGModel;
    if (trgVTK) delete trgVTK;
  }
  else
  {  
    std::cout << "WORKING WITH DIM=1" << std::endl;
    adapter->registerData(array_name, flatten(pntData));
    // check if data is correctly attached
    std::vector<double> preAMRData;
    adapter->getMeshData(array_name, &preAMRData);
    std::cout << "preAMRData[1]  = " << preAMRData[1] << std::endl;
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
    MAd::M_writeMsh(mesh, (char*) &outMeshFile[0u], 2);
   
    // get data after refinement
    std::vector<double> postAMRData;
    adapter->getMeshData(array_name, &postAMRData);
    std::cout << "postAMRData[1] = " << postAMRData[1] << std::endl;
    std::cout << "Size of vector = " << postAMRData.size() << std::endl; 

    // convert refined msh file back to vtk
    GModel* trgGModel;
    trgGModel = new GModel("refined"); 
    trgGModel->readMSH((char*) &outMeshFile[0u]);
    trgGModel->writeVTK(trim_fname(outMeshFile,".vtk"), false, true); // binary=false, saveall=true


    // write physical quantities to vtk file
    vtkAnalyzer* trgVTK;
    trgVTK = new vtkAnalyzer((char*) &(trim_fname(outMeshFile,".vtk")[0u]));
    trgVTK->read();
    trgVTK->setPointDataArray(&array_name[0u], 1, postAMRData);
    trgVTK->report();
    trgVTK->write((char*) &(trim_fname(outMeshFile, "_solution.vtu"))[0u]);
    trgVTK->writeMSH(trim_fname(outMeshFile, "_solution.msh"));
    if (trgGModel) delete trgGModel;
    if (trgVTK) delete trgVTK;
  }

 
  //if (localSF) delete localSF;
  //if (adapter) delete adapter;
  if (bSF) delete bSF;

  return 0;
}

//-------------------------- Auxiliary Functions --------------------------------//

// computes reciprocal of argument
template <typename T>
T reciprocal(T x)
{
  return (T) 1/x ;
}

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x)
{
  std::transform(x.begin(), x.end(), x.begin(), reciprocal<T>);
}

// find minmax of vector, excluding inf 
std::vector<double> getMinMax(const std::vector<double>& x)
{
  std::vector<double> result(2);
  std::vector<double> tmp;
  for (int i = 0; i < x.size(); ++i)
  {
    if (!std::isinf(x[i])) // exclude inf
    {
      tmp.push_back(x[i]); 
    }
  }
  auto minmax = std::minmax_element(tmp.begin(), tmp.end());
  result[0] = *minmax.first;
  result[1] = *minmax.second;
  return result;
} 

// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
double scale_to_range(double x, const std::vector<double>& xminmax, 
                                const std::vector<double>& yminmax)
{
  if (std::isinf(x))
    return yminmax[1];
  return yminmax[0] + (yminmax[1] - yminmax[0])*(x - xminmax[0])/(xminmax[1] - xminmax[0]);
}

// scales each number in vector in range [xmin, xmax] to [ymin,ymax]
void scale_vec_to_range(std::vector<double>& x, 
                        const std::vector<double>& xminmax,
                        const std::vector<double>& yminmax)
{
  for (int i = 0; i < x.size(); ++i)
    x[i] = scale_to_range(x[i], xminmax, yminmax);
}

// multiplies vector by scalar, in place
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x)
{
  std::vector<T> result;
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), result.begin(), std::bind1st(std::multiplies<T>(),a));
  return result;
}  

// get average and stdev of values
std::vector<double> getMeanStdev(std::vector<double>& x)
{
  std::vector<double> result(2);
  double ave = 0;
  for (int i = 0; i < x.size(); ++i)
    ave += x[i];
  ave /= x.size();

  double stdev = 0;
  for (int i = 0; i < x.size(); ++i)
    stdev += (x[i] - ave)*(x[i] - ave);
  stdev = std::sqrt(stdev/x.size());
  result[0] = ave;
  result[1] = stdev;
  return result;
}

// generates boolean array with 1 if value >= tol, 0 otherwise
std::vector<int> cellsToRefine(std::vector<double>& values, double tol)
{
  std::vector<int> result(values.size(),0);  
  for (int i = 0; i < values.size(); ++i)
  {
    if (values[i] > tol)
      result[i] = 1;
  }
  return result;
}

std::string trim_fname(std::string fname, std::string ext)
{
  size_t beg = 0;
  size_t end = fname.find(".");
  std::string name;
  if (end != -1)
  {
    name = fname.substr(beg,end);
    name.append(ext);
    return name;
  }

  else 
  {
    std::cout << "Error finding file extension for " << fname << std::endl;
    exit(1);
  }
}  
//---------------------------------------------------------------------------//

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
