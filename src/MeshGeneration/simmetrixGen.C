#include <iostream>
#include <simmetrixGen.H>
#include <simmetrixParams.H>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>


// constructor
simmetrixGen::simmetrixGen(simmetrixParams* _params)
    : params(_params), haveLog(false), writeSurfAndVol(false),
      prog(NULL), model(NULL), dModel(NULL), mcase(NULL), mesh(NULL)
{

  if (params->logFName != "NONE")
  {
      haveLog = true;
      Sim_logOn(&(params->logFName)[0u]);
  }

  std::cout << "Log File: " << params->logFName << std::endl;  
  std::cout << "License File: " << params->licFName << std::endl;  
  std::cout << "Feature List: " << params->features << std::endl;  

  // initialization
  SimLicense_start(&(params->features)[0u], &(params->licFName)[0u]);
  MS_init();
  Sim_setMessageHandler(messageHandler);
  setProgress();
  std::cout << "Simmetrix mesh generator created" << std::endl;
}

// default
simmetrixGen::simmetrixGen()
    : params(NULL), haveLog(true), writeSurfAndVol(false),
      prog(NULL), model(NULL), dModel(NULL), mcase(NULL), mesh(NULL)
{
  Sim_logOn("simmetrixGen.log");
  SimLicense_start("geomsim_core,meshsim_surface,meshsim_volume","simmodsuite.lic");
  MS_init();
  Sim_setMessageHandler(messageHandler);
  setProgress(); 
  std::cout << "Simmetrix mesh generator created" << std::endl;
}

// destructor
simmetrixGen::~simmetrixGen()
{
  if (mesh)
    M_release(mesh);
  if (mcase)
    MS_deleteMeshCase(mcase);
  if (model)
    GM_release(model);
  if (dModel)
    GM_release(dModel);
  if (prog)
    Progress_delete(prog);    
  SimLicense_stop();
  MS_exit();
  if (haveLog)
    Sim_logOff();
  std::cout << "Simmetrix mesh generator destroyed" << std::endl;
}

void simmetrixGen::setProgress()
{
  prog = Progress_new();
  Progress_setDefaultCallback(prog);    
}

int simmetrixGen::createSurfaceMeshFromSTL(const char* stlFName)
{
  if (!createModelFromSTL(stlFName))
  {
    try
    {
      SimDiscrete_start(0);
      pModelItem modelDomain = GM_domain(dModel);
      mcase = MS_newMeshCase(dModel);
      std::cout << "Number of initial mesh faces: " << M_numFaces(mesh) << std::endl;

      if (MS_checkMeshIntersections(mesh,0,prog))
      {
        std::cerr << "There are intersections in the input surface mesh" << std::endl;
        MS_deleteMeshCase(mcase);
        M_release(mesh);
        GM_release(dModel);
        return 1;
      }

      // Mesh Size
      //--------------------------------------------
      // set the mesh size
      // params->meshSize is multiplied by the largest edge of the coordinate-aligned bonding box
      // to give a typical size of the cell edges
      MS_setMeshSize(mcase, modelDomain, 2, params->meshSize, 0);

      // Anisotropic Curvature Refinement
      //--------------------------------------------
      // setting curvature-based refinement params
      // params->anisoMeshCurv must be < 0.5
      // The mesh size h is selected such that d/h < par. par should always be less than 0.5.
      // Useful values for par are typically in the range of 0.01 to 0.4 (smaller value = more refinement).
      MS_setAnisoMeshCurv(mcase, modelDomain, 2, params->anisoMeshCurv);

      // Minimum Curvature Size (not set in JSON, defaults to 0.01)
      //--------------------------------------------
      // params->minCurvSize must be < 0.5
      // The mesh size h is selected such that d/h < par. par should always be less than 0.5.
      // Useful values for par are typically in the range of 0.01 to 0.4 (smaller value = more refinement).
      MS_setMinCurvSize(mcase, modelDomain, 2, params->minCurvSize);

      // Global Gradation Rate
      //--------------------------------------------
      // Allows faster transition of mesh size from finer to coarser areas.
      // By default, mesh size approximately doubles with each element until it reaches the coarser size.
      // This function allows a faster transition. For a slower than default transition, MS_setGlobalSizeGradationRate() should be used.
      // Valid values of factor are between 0 and 1. A factor of 1 is the same as default, ie. the mesh size will double with each element.
      // As factor decreases, the transition in size will happen in fewer elements (larger jumps in element size), until at a factor = 0,
      // the mesh will immediately jump to the coarser mesh size with no transition.
      MS_setGlobalSizeGradationRate(mcase, params->glbSizeGradRate);

      pSurfaceMesher surfaceMesher = SurfaceMesher_new(mcase, mesh);
      // snap model vertices to resulting surface mesh (ensures consistency with model)
      // this causes seg faults when trying to change topology of mesh during improvement
      //SurfaceMesher_setSnapForDiscrete(surfaceMesher, 1);
      SurfaceMesher_setEnforceSpatialGradation(surfaceMesher, 1);
      SurfaceMesher_execute(surfaceMesher, prog);
      std::cout << "Number of mesh faces on the surface: " << M_numFaces(mesh) << std::endl;
      SurfaceMesher_delete(surfaceMesher);
      // improving the surface mesh
      pSurfaceMeshImprover surfaceMeshImprover = SurfaceMeshImprover_new(mesh);
      SurfaceMeshImprover_setShapeMetric(surfaceMeshImprover, ShapeMetricType_MaxAngle,145.);
      // set smoothing type to move vertices down gradient of shape metric
      SurfaceMeshImprover_setSmoothType(surfaceMeshImprover, 1);

      // Surface Mesh Improver Gradation Rate
      //--------------------------------------------
      // Sets the mesh size gradation rate for improvement refinement.
      // In order to meet the shape metric criteria, the surface mesh improver may locally refine the mesh.
      // This control sets the gradation rate from the refined entities. Similar to setGlobalSizeGradationRate.
      // Default value is 2/3.
      SurfaceMeshImprover_setGradationRate(surfaceMeshImprover, params->surfMshImprovGradRate);


      // Surface Mesh Improver Min Size
      //--------------------------------------------
      // Allows the mesh improver to refine the mesh to meet the specified criteria.
      // If this function is called, the mesh improver will be permitted to refine the mesh in areas where
      // required to meet the specified criteria. The size passed in is the minimum size to which the mesh
      // will be refined.
      SurfaceMeshImprover_setMinRefinementSize
        (surfaceMeshImprover, 2, params->surfMshImprovMinSize);

      // fix intersections after improvement
      SurfaceMeshImprover_setFixIntersections(surfaceMeshImprover, 2); 
      SurfaceMeshImprover_execute(surfaceMeshImprover, prog);
      SurfaceMeshImprover_delete(surfaceMeshImprover);
      std::cout << "Num regions in surf mesh: " << M_numRegions(mesh) << std::endl;
      SimDiscrete_stop(0);
    }
    catch (pSimError err) 
    {
      std::cerr << "SimModSuite error caught:" << std::endl;
      std::cerr << "  Error code: " << SimError_code(err) << std::endl;
      std::cerr << "  Error string: " << SimError_toString(err) << std::endl;
      SimError_delete(err);
      return 1;
    } 
    catch (...) 
    {
      std::cerr << "Unhandled exception caught during surface mesh generation" << std::endl;
      return 1;
    }
    return 0; 
  }
  else
    return 1;
}

int simmetrixGen::createVolumeMeshFromSTL(const char* stlFName)
{
  if (!createSurfaceMeshFromSTL(stlFName))
  {
    try
    {
      SimDiscrete_start(0);
      pVolumeMesher volumeMesher = VolumeMesher_new(mcase, mesh);
      VolumeMesher_setEnforceSize(volumeMesher, 1);
      VolumeMesher_execute(volumeMesher, prog);
      VolumeMesher_delete(volumeMesher);
      SimDiscrete_stop(0);
    }
    catch (pSimError err) 
    {
      std::cerr << "SimModSuite error caught:" << std::endl;
      std::cerr << "  Error code: " << SimError_code(err) << std::endl;
      std::cerr << "  Error string: " << SimError_toString(err) << std::endl;
      SimError_delete(err);
      return 1;
    } 
    catch (...) 
    {
      std::cerr << "Unhandled exception caught during volume mesh generation" << std::endl;
      return 1;
    }
    return 0; 
  }
  else
    return 1;
}

int simmetrixGen::createMeshFromSTL(const char* stlFName)
{
  if (!createVolumeMeshFromSTL(stlFName))
  {
    convertToVTU();
    return 0;
  }
  return 1;
}

int simmetrixGen::createModelFromSTL(const char* stlFName)
{
  // try/catch around all SimModSuite calls
  // as errors are thrown.
  try 
  { 

    SimDiscrete_start(0); 
    mesh = M_new(0,0);
    //pDiscreteModel model = 0;
    // load and return if error encountered
    if(M_importFromSTLFile(mesh, stlFName, prog)) 
    { 
      std::cerr << "Error importing file" << std::endl;
      M_release(mesh);
      return 1;
    }
  
    // check the input mesh for intersections
    // this call must occur before the discrete model is created
    if(MS_checkMeshIntersections(mesh,0,prog)) 
    {
      std::cerr << "There are intersections in the input mesh" << std::endl;
      M_release(mesh);
      return 1;
    }
    
    // create the Discrete model
    dModel = DM_createFromMesh(mesh, 1, prog);
    // return if error in model creation
    if(!dModel) 
    { 
      std::cerr << "Error creating Discrete model from mesh" << std::endl;
      M_release(mesh);
      return 1;
    }
    
    // define the Discrete model
    DM_findEdgesByFaceNormals(dModel, 5., prog);
    DM_eliminateDanglingEdges(dModel, prog);
    // complete the topology and return if erro encountered
    if(DM_completeTopology(dModel, prog)) 
    { 
      std::cerr << "Error completing Discrete model topology" << std::endl;
      M_release(mesh);
      GM_release(dModel);
      return 1;
    }
    
    // Print out information about the model
    std::cout << "Number of model vertices: " << GM_numVertices(dModel) << std::endl;
    std::cout << "Number of model edges: " << GM_numEdges(dModel) << std::endl;
    std::cout << "Number of model faces: " << GM_numFaces(dModel) << std::endl;
    std::cout << "Number of model regions: " << GM_numRegions(dModel) << std::endl;
    
    GM_write(dModel,"discreteModel.smd",0,prog); // save the discrete model
    SimDiscrete_stop(0);
  } 
  catch (pSimError err) 
  {
    std::cerr << "SimModSuite error caught:" << std::endl;
    std::cerr << "  Error code: " << SimError_code(err) << std::endl;
    std::cerr << "  Error string: " << SimError_toString(err) << std::endl;
    SimError_delete(err);
    return 1;
  } 
  catch (...) 
  {
    std::cerr << "Unhandled exception caught during discrete model generation step" << std::endl;
    return 1;
  }
  return 0; 
}

void simmetrixGen::convertToVTU()
{
  if (!dataSet)
  {
    // points to be pushed into dataSet
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // declare vtk dataset
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp 
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    // get the iterator for mesh vertices 
    VIter vertices = M_vertexIter(mesh);
    pVertex vertex;
    double coord[3];
    // allocate size of vtk point container
    points->SetNumberOfPoints(M_numVertices(mesh));
    int i = 0;
    // copy points into vtk point container
    vertex = VIter_next(vertices);
    while (vertex)
    {
      EN_setID( (pEntity) vertex, i); 
      V_coord(vertex, coord);
      points->SetPoint(i,coord);
      ++i;
      vertex = VIter_next(vertices);
    }
    VIter_delete(vertices);
    dataSet_tmp->SetPoints(points);

    // if mesh is a surface mesh
    if (!M_numRegions(mesh))
    {
      // allocate space for faces
      dataSet_tmp->Allocate(M_numFaces(mesh));
      FIter faces = M_faceIter(mesh);
      pFace face;
      pPList faceVerts;
      face = FIter_next(faces);
      while (face)
      {
        faceVerts = F_vertices(face,1);
        createVtkCell(dataSet_tmp, 3, VTK_TRIANGLE, faceVerts);
        PList_delete(faceVerts);
        face = FIter_next(faces);
      }
      FIter_delete(faces);
      dataSet = dataSet_tmp;
    }
    // if mesh is volume mesh and you want to write surf and vol cells
    else if (writeSurfAndVol)
    {
      // allocate space for cells in vtk dataset
      dataSet_tmp->Allocate(M_numRegions(mesh) + M_numFaces(mesh));
      // add surface faces
      FIter faces = M_faceIter(mesh);
      pFace face;
      pPList faceVerts;
      face = FIter_next(faces);
      while (face)
      {
        faceVerts = F_vertices(face,1);
        createVtkCell(dataSet_tmp, 3, VTK_TRIANGLE, faceVerts);
        PList_delete(faceVerts);
        face = FIter_next(faces);
      }
      FIter_delete(faces);
      // add vol cells
      addVtkVolCells(dataSet_tmp);
      dataSet = dataSet_tmp;
    }
    // write only vol cells
    else
    {
      dataSet_tmp->Allocate(M_numRegions(mesh)); 
      addVtkVolCells(dataSet_tmp);
      dataSet = dataSet_tmp;
    }
  }
}

void simmetrixGen::setWriteSurfAndVol(bool b)
{
  writeSurfAndVol = b;
}

void simmetrixGen::saveMesh(const std::string& mFName)
{

  size_t last = mFName.find_last_of('.');
  if (last != std::string::npos)
  {
    std::string ext = mFName.substr(last); 
    if (ext == ".sms")
    {
      std::cout << "writing .sms mesh" << std::endl;
      M_write(mesh, &mFName[0u],0,prog);
    } 
    else if (ext == ".vtu")
    {
      std::cout << "writing .vtu mesh" << std::endl;
      convertToVTU();
      writeVTFile<vtkXMLUnstructuredGridWriter>(mFName,dataSet);  
    } 
  }
}

void simmetrixGen::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                              const int numIds,
                              const int cellType,
                              pPList regionVerts)
{
  pVertex vertex;
  vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkCellIds->SetNumberOfIds(numIds);
  for (int i = 0; i < numIds; ++i)
  {
    vertex = (pVertex) PList_item(regionVerts,i);
    vtkCellIds->SetId(i,EN_id( (pEntity) vertex));
  }
  dataSet->InsertNextCell(cellType,vtkCellIds);
}
                             
void simmetrixGen::addVtkVolCells(vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp)
{
  // get iterator for mesh regions (cells) and add them
  RIter regions = M_regionIter(mesh);
  pRegion region;  
  pPList regionVerts;
  region = RIter_next(regions);
  while (region)
  {
    regionVerts = R_vertices(region,0);
    rType celltype = R_topoType(region);
    switch(celltype)
    {
      case Rtet:
      {
        createVtkCell(dataSet_tmp, 4, VTK_TETRA, regionVerts);
        break;
      }
      case Rwedge:
      {
        createVtkCell(dataSet_tmp, 5, VTK_WEDGE, regionVerts);
        break;
      }
      case Rpyramid:
      {
        createVtkCell(dataSet_tmp, 4, VTK_PYRAMID, regionVerts);
        break;
      }
      case Rhex:
      {
        createVtkCell(dataSet_tmp, 6, VTK_HEXAHEDRON, regionVerts);
        break;
      }
      case Rdwedge:
      case Rdhex:
      case Runknown:
        std::cerr << "Encountered unknown cell type: " << celltype << std::endl;
        exit(1);
    }
    PList_delete(regionVerts);
    region = RIter_next(regions);
  }
  RIter_delete(regions);
}

void simmetrixGen::messageHandler(int type, const char* msg)
{
  switch (type) 
  {
    case Sim_InfoMsg:
      std::cout << "Info: " << msg <<std::endl;
      break;
    case Sim_DebugMsg:
      std::cout << "Debug: " << msg <<std::endl;
      break;
    case Sim_WarningMsg:
      std::cout << "Warning: " << msg << std::endl;
      break;
    case Sim_ErrorMsg:
      std::cout << "Error: " << msg<< std::endl;
      break;
  }
}

void simmetrixGen::createMeshFromModel(const char* mFName)
{
  model = GM_load(mFName, 0, prog);
  mcase = MS_newMeshCase(model);
  MS_setMeshSize(mcase, GM_domain(model), 2, 0.5, 0);
  mesh = M_new(0, model);
  pSurfaceMesher surfaceMesher = SurfaceMesher_new(mcase, mesh);
  SurfaceMesher_execute(surfaceMesher, prog);
  SurfaceMesher_delete(surfaceMesher);
  pVolumeMesher volumeMesher = VolumeMesher_new(mcase, mesh);
  VolumeMesher_execute(volumeMesher, prog);
  VolumeMesher_delete(volumeMesher);
}
