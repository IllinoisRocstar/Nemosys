// Nemosys
#include <meshPhys.H>
#include <netgenInterface.H>
// stl
#include <chrono>

/******************************************************************************
* Usage: refineMesh meshFile outputMesh array_id refine_method stdev_mult     *
                                                                              * 
* meshFile is the mesh to be refined                                          *
* outputMesh is the .msh name of the refined mesh                             *
* array_id is the data on which to base the refinement                        *
   - use `xmlDump meshFile` for list of available ids                         *
* refine_method is "val" (value based) or "grad" (gradient based)             *
* stdev_mult is the multiplier on stdev of data that defines                  *
  a threshold value. Cells with data values larger than this will be refined  *
*******************************************************************************

TODO:  writing vector data to pos file not implemented in madlib.
       instead, split vector data into component vectors
       if writing data required. This isn't required for anything though.
TODO:  Register all data with the nodal data manager so that output refined
       mesh contains everything 
TODO:  GModel::writeVTK method only writes points and connectivities to vtk,
       no data. Possible class hierarchy could begin with providing interface 
       between all formats (msh, vtk, netgen) for 
          1) points
          2) cells and correct connectivities (differ in order by format)
          3) point data
          4) cell data
          5) metadata like material tags etc.
       Then, we store all of these in one Mesh class and apply various methods
       based on what is available with the format. 
******************************************************************************/

//---------------------------Auxiliary Classes---------------------------------//

// class wrapping around chrono for timing methods
class Timer {
private:
  typedef std::chrono::time_point<std::chrono::system_clock> time_t;

public:
  Timer() : startTime(), stopTime() {}

  time_t start()   { return (startTime = std::chrono::system_clock::now()); }
  time_t stop()    { return (stopTime  = std::chrono::system_clock::now()); }
  double elapsed() { return std::chrono::duration_cast<std::chrono::milliseconds>
                                                      (stopTime-startTime).count(); }

private:
  time_t startTime, stopTime;
};

//----------------------------------------------------------------------//

int main(int argc, char* argv[])
{
//  // Check input
//  if (argc < 7  && strcmp(argv[1], "Nek")) 
//  {
//    std::cout << "Usage Error: "  << argv[0] 
//                                  << " meshFile outputMesh array_id"
//                                  << " refine_method stdev_mult maxIsmin <write_test>\n \n" 
//                                  << "where <write_test> is an optional bool (1|0) for testing"
//                                  << std::endl;
//    exit(1);
//  }
//
//  // seperate orchestration for nek stuff
//  if (!strcmp(argv[1], "Nek"))
//  {
//    // input mesh
//    meshPhys* mshphys;
//    mshphys = new meshPhys(argv[2]);
//    mshphys->report();
//    // method to determine cells to refine ("grad", "val")
//    std::string method = argv[3];
//    // multiplier for stdev of data to filter cells considered for refinement
//    double dev_mult = std::stof(argv[4]);
//    std::cout << "Writing cells to refine for Nek ..." << std::endl;
//    mshphys->writeCellsToRefine(1, method, dev_mult);
//    std::cout << "Success!" << std::endl; 
//    delete mshphys;
//    return 0;
//  }
//  
//  // input mesh
//  std::string meshFile = argv[1];
//  // output mesh
//  std::string outMeshFile = argv[2];
//  // data array index
//  int array_id = std::stoi(argv[3]);
//  // method to generate sizes ("grad", "val")
//  std::string method = argv[4];
//  // multiplier for stdev of data to filter cells considered for refinement
//  double dev_mult = std::stof(argv[5]);
//  // boolean for whether max elem length after refinement is min or multiplier on max
//  bool maxIsmin = (bool) std::stoi(argv[6]);
// 
//  // create gmodel from input mesh
//  MAd::pGModel gmodel = 0;
//  MAd::pMesh mesh = MAd::M_new(gmodel);
//  meshPhys* mshphys;
//  
//  // check for vtk extension
//  if (meshFile.find(".v") != -1 ) 
//  {
//    mshphys = new meshPhys((char*) &meshFile[0u]);
//    std::string ofname = "converted.msh";
//    // convert vtk to msh for loading into MAdLib format 
//    mshphys->writeMSH(ofname);  
//    MAd::M_load(mesh, "converted.msh"); 
//  }
//  // check for gmsh extension
//  else if (meshFile.find(".msh") != -1) 
//  {
//    MAd::M_load(mesh, (char*) &(meshFile)[0u]); 
//    GModel* tmpGModel;
//    tmpGModel = new GModel("tmp");
//    tmpGModel->readMSH((char*) &meshFile[0u]);
//    std::string new_name = trim_fname(meshFile,".vtk");
//    // convert msh to vtk for instantiation of meshPhys
//    tmpGModel->writeVTK(new_name, false, true); // binary=false, saveall=true
//    delete tmpGModel;
//    mshphys = new meshPhys((char*) &new_name[0u]);
//  }
//
//  // define a size field over the mesh and write to backgroundSF.msh 
//  mshphys->createSizeField(array_id, method, dev_mult, maxIsmin); 
//  // extracting dimension (values per node) for considered array
//  int dim = mshphys->getDimArray(array_id);
//  // print a report of the vtk mesh
//  mshphys->report();
//
//  // finding and classifying boundary elements as 2d 
//  mesh->classify_unclassified_entities();
//  mesh->destroyStandAloneEntities();
//  MAd::pGEntity bnd = (MAd::pGEntity) MAd::GM_faceByTag(mesh->model, 0);
//  mesh->classify_grid_boundaries(bnd);
//
//  // loading size field into background sizefield instantiation
//  MAd::BackgroundSF* bSF = new MAd::BackgroundSF("backgroundSF");
//  bSF->loadData("backgroundSF.msh");
//
//  std::cout << "\n \n Beginning Adapter Construction" << std::endl;
//  // timing adapter construction
//  Timer T;
//  T.start();
//  // instantiating adapter with background sizefield
//  MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, bSF);
//  T.stop();
//  std::cout << "Time for adapter construction (ms): " << T.elapsed() << "\n \n";
// 
//  // RUNNING ADAPTER TO REFINE MESH 
//  mshphys->Refine(adapter, mesh, array_id, dim, outMeshFile);    
//  /* This function essentially does the following:
//      1) Registers the relevant data with the nodal data manager
//      2) Runs the adapter
//      3) Unclassifies the boundary elements for proper output
//      4) Writes the refined mesh to a .msh file without data
//      5) Gets the data after refinement
//      6) Converts the refined mesh to a vtk file without data
//      7) Writes the post-refinement data to a refined vtk and gmsh mesh 
//      8) Misc: Memory management, string trimming/file naming      */
// 
//  //--------------- write test files if testing enabled --------------------------------//
//  bool write_test = 0;
//  if (argc == 8)
//    write_test = (bool) std::stoi(argv[7]);
//  
//  if (write_test) 
//  {
//    // extracting array_name for proper naming in output
//    std::string array_name =  mshphys->getPointData(array_id).getName(); 
//    if (method.compare("grad") == 0)
//    {
//      std::vector<double> result = mshphys->ComputeL2GradAtAllCells(array_id);
//      array_name.insert(0,"L2Grad_");
//      mshphys->setCellDataArray(&array_name[0u], 1, result); 
//      array_name.append("AtCell.vtu");
//      mshphys->write((char*) &array_name[0u]);
//      std::ofstream outputStream("GradTest.txt");
//      if (!outputStream.good())
//      {
//        std::cout << "error opening GradTest.txt" << std::endl;
//        exit(2);
//      }
//      for (int i = 0; i < result.size(); ++i)
//        outputStream << result[i] << std::endl;
//    }
//    else if (method.compare("val") == 0)
//    {
//      std::vector<double> result = mshphys->ComputeL2ValAtAllCells(array_id);
//      array_name.insert(0,"L2Val_");
//      mshphys->setCellDataArray(&array_name[0u], 1, result); 
//      array_name.append("AtCell.vtu");
//      mshphys->write((char*) &array_name[0u]);
//      std::ofstream outputStream("ValTest.txt");
//      if (!outputStream.good())
//      {
//        std::cout << "error opening GradTest.txt" << std::endl;
//        exit(2);
//      }
//      for (int i = 0; i < result.size(); ++i)
//        outputStream << result[i] << std::endl;
//    }
//    else
//    {
//      std::cout << "method must be \"val\" or \"grad\" " << std::endl;
//      exit(2);
//    }
//  }
//
//  if (bSF) delete bSF;
//  if (mshphys) delete mshphys;
  netgenInterface* tmp = new netgenInterface();
  tmp->createMeshFromSTL("hinge.stl");   
  tmp->exportToVTK("hinge.vtk");
  netgenInterface* tmp1 = new netgenInterface(); 
  tmp1->importFromVol("hinge.vol");
  tmp1->exportToVTK("hinge1.vtk");
  //Ng_SolutionData* tmp1;
  //Ng_InitSolutionData (tmp1);
 // tmp->importFromVTK(argv[1]);
  
  if(tmp) delete tmp;
  if(tmp1) delete tmp1;
  return 0;
}

//-----------------------------------------------------------------------------------------//

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
  /*mshphys->setCellDataArray(&array_name[0u], dim, values); 
  mshphys->report();
  std::string newname;
  newname.insert(0,array_name);
  newname.insert(0, "_");
  newname.append("AtCell.vtu");
  mshphys->write((char*) &(trim_fname(meshFile, newname))[0u]);
 */
