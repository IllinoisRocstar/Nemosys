/* 
   The objective of this utility is to test Nemosys interface to DTK data transfer
   tools.
*/


// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include <cgnsAnalyzer.H>
#include <vtkAnalyzer.H>
#include <baseInterp.H>
#include <meshPartitioner.H>
#include <cgnsWriter.H>

// MAdLib headers 
#include <ModelInterface.h>
#include <MAdLib.h>
#include <NodalDataManager.h>
#include <MeshDataBaseIO.h>

// Others
#include "DTK_STKMeshHelpers.hpp"
#include "DTK_STKMeshManager.hpp"
#include "DTK_MapOperatorFactory.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Tpetra_MultiVector.hpp>

#include <Intrepid_FieldContainer.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_topology/topology.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>



// typedefs
//typedef std::shared_ptr<cgnsAnalyzer> cgPtr;


/* auxiliary functions */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds);
std::vector<double> getCOG(std::vector<double>& regCntCrds);
double dataFunction( double x, double y, double z );

/*   Main Function */ 
int main(int argc, char* argv[])
{
    // INITIALIZATION
    // --------------

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
   
    // Read in command line options.
    std::string xml_input_filename;
    Teuchos::CommandLineProcessor clp(false);
    clp.setOption( "xml-in-file",
                   &xml_input_filename,
                   "The XML file to read into a parameter list" );
    clp.parse(argc,argv);

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
        Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
        xml_input_filename, Teuchos::inoutArg(*plist) );

    // Read command-line options
    std::string source_mesh_input_file =
        plist->get<std::string>("Source Mesh Input File");
    std::string source_mesh_output_file =
        plist->get<std::string>("Source Mesh Output File");
    std::string source_mesh_part_name =
        plist->get<std::string>("Source Mesh Part");
    std::string target_mesh_input_file =
        plist->get<std::string>("Target Mesh Input File");
    std::string target_mesh_output_file =
        plist->get<std::string>("Target Mesh Output File");
    std::string target_mesh_part_name =
        plist->get<std::string>("Target Mesh Part");

    // Get the raw mpi communicator (basic typedef in STK).
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm =
        mpi_comm->getRawMpiComm();
    stk::ParallelMachine parallel_machine = (*opaque_comm)();


    // SOURCE MESH READ
    // ----------------

    // Load the source mesh.
    stk::io::StkMeshIoBroker src_broker( parallel_machine );
    std::size_t src_input_index = src_broker.add_mesh_database(
        source_mesh_input_file, "exodus", stk::io::READ_MESH );
    src_broker.set_active_mesh( src_input_index );
    src_broker.create_input_mesh();

    // Add a nodal field to the source part.
    //std::cout << src_broker.meta_data().get_part(0).name() << std::endl;
    stk::mesh::Field<double>& source_field =
        src_broker.meta_data().declare_field<stk::mesh::Field<double> >(
            stk::topology::NODE_RANK, "u_src" );
    stk::mesh::Part* src_part =
        src_broker.meta_data().get_part( source_mesh_part_name );
    if (src_part == NULL){
      std::cout << "Source part not found\n";
      throw;
    }
    
    stk::mesh::put_field( source_field, *src_part );

    // Create the source bulk data.
    src_broker.populate_bulk_data();
    Teuchos::RCP<stk::mesh::BulkData> src_bulk_data =
        Teuchos::rcpFromRef( src_broker.bulk_data() );

    // Put some data in the source field. We will use the distance of the node
    // from the origin as the data.
    stk::mesh::Selector src_stk_selector( *src_part );
    stk::mesh::BucketVector src_part_buckets =
        src_stk_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> src_part_nodes;
    stk::mesh::get_selected_entities(
        src_stk_selector, src_part_buckets, src_part_nodes );
    Intrepid::FieldContainer<double> src_node_coords =
        DataTransferKit::STKMeshHelpers::getEntityNodeCoordinates(
            Teuchos::Array<stk::mesh::Entity>(src_part_nodes), *src_bulk_data );
    double* src_field_data;
    int num_src_part_nodes = src_part_nodes.size();
    for ( int n = 0; n < num_src_part_nodes; ++n )
    {
        src_field_data = stk::mesh::field_data( source_field, src_part_nodes[n] );
        src_field_data[0] = dataFunction( src_node_coords(n,0,0),
                                          src_node_coords(n,0,1),
                                          src_node_coords(n,0,2) );
    }

    // TARGET MESH READ
    // ----------------

    // Load the target mesh.
    stk::io::StkMeshIoBroker tgt_broker( parallel_machine );
    std::size_t tgt_input_index = tgt_broker.add_mesh_database(
        target_mesh_input_file, "exodus", stk::io::READ_MESH );
    tgt_broker.set_active_mesh( tgt_input_index );
    tgt_broker.create_input_mesh();

    // Add a nodal field to the target part.
    stk::mesh::Field<double>& target_field =
        tgt_broker.meta_data().declare_field<stk::mesh::Field<double> >(
            stk::topology::NODE_RANK, "u_tgt" );
    stk::mesh::Part* tgt_part =
        tgt_broker.meta_data().get_part( target_mesh_part_name );
    stk::mesh::put_field( target_field, *tgt_part );

    // Add an error nodal field to the target part.
    stk::mesh::Field<double>& target_error_field =
        tgt_broker.meta_data().declare_field<stk::mesh::Field<double> >(
            stk::topology::NODE_RANK, "u_err" );
    stk::mesh::put_field( target_error_field, *tgt_part );

    // Create the target bulk data.
    tgt_broker.populate_bulk_data();
    Teuchos::RCP<stk::mesh::BulkData> tgt_bulk_data =
        Teuchos::rcpFromRef( tgt_broker.bulk_data() );

    // SOLUTION TRANSFER SETUP
    // -----------------------

    // Create a manager for the source part elements.
    DataTransferKit::STKMeshManager src_manager( src_bulk_data, src_stk_selector );

    // Create a manager for the target part nodes.
    stk::mesh::Selector tgt_stk_selector( *tgt_part );
    DataTransferKit::STKMeshManager tgt_manager( tgt_bulk_data, tgt_stk_selector );

    // Create a solution vector for the source.
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > src_vector =
        src_manager.createFieldMultiVector<stk::mesh::Field<double> >(
            Teuchos::ptr(&source_field), 1 );

    // Create a solution vector for the target.
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > tgt_vector =
        tgt_manager.createFieldMultiVector<stk::mesh::Field<double> >(
            Teuchos::ptr(&target_field), 1 );

    // Print out source mesh info.
    Teuchos::RCP<Teuchos::Describable> src_describe =
        src_manager.functionSpace()->entitySet();
    std::cout << "Source Mesh" << std::endl;
    src_describe->describe( std::cout );
    std::cout << std::endl;

    // Print out target mesh info.
    Teuchos::RCP<Teuchos::Describable> tgt_describe =
        tgt_manager.functionSpace()->entitySet();
    std::cout << "Target Mesh" << std::endl;
    tgt_describe->describe( std::cout );
    std::cout << std::endl;

    // SOLUTION TRANSFER
    // -----------------

    // Create a map operator. The operator settings are in the
    // "DataTransferKit" parameter list.
    Teuchos::ParameterList& dtk_list = plist->sublist("DataTransferKit");
    DataTransferKit::MapOperatorFactory op_factory;
    Teuchos::RCP<DataTransferKit::MapOperator> map_op =
        op_factory.create( src_vector->getMap(),
                           tgt_vector->getMap(),
                           dtk_list );

    // Setup the map operator. This creates the underlying linear operators.
    map_op->setup( src_manager.functionSpace(), tgt_manager.functionSpace() );

    // Apply the map operator. This interpolates the data from one STK field
    // to the other.
    map_op->apply( *src_vector, *tgt_vector );

    // COMPUTE THE SOLUTION ERROR
    // --------------------------

    stk::mesh::BucketVector tgt_part_buckets =
        tgt_stk_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> tgt_part_nodes;
    stk::mesh::get_selected_entities(
        tgt_stk_selector, tgt_part_buckets, tgt_part_nodes );
    Intrepid::FieldContainer<double> tgt_node_coords =
        DataTransferKit::STKMeshHelpers::getEntityNodeCoordinates(
            Teuchos::Array<stk::mesh::Entity>(tgt_part_nodes), *tgt_bulk_data );
    double* tgt_field_data;
    double* err_field_data;
    int num_tgt_part_nodes = tgt_part_nodes.size();
    double error_l2_norm = 0.0;
    double field_l2_norm = 0.0;
    for ( int n = 0; n < num_tgt_part_nodes; ++n )
    {
        double gold_value = dataFunction( tgt_node_coords(n,0,0),
                                          tgt_node_coords(n,0,1),
                                          tgt_node_coords(n,0,2) );
        tgt_field_data = stk::mesh::field_data( target_field, tgt_part_nodes[n] );
        err_field_data = stk::mesh::field_data( target_error_field, tgt_part_nodes[n] );
        err_field_data[0] = tgt_field_data[0] - gold_value;
        error_l2_norm += err_field_data[0] * err_field_data[0];
        field_l2_norm += tgt_field_data[0] * tgt_field_data[0];
        err_field_data[0] /= gold_value;
    }
    error_l2_norm = std::sqrt( error_l2_norm );
    field_l2_norm = std::sqrt( field_l2_norm );
    std::cout << "|e|_2 / |f|_2: " << error_l2_norm / field_l2_norm << std::endl;


    // SOURCE MESH WRITE
    // -----------------

    std::size_t src_output_index = src_broker.create_output_mesh(
        source_mesh_output_file, stk::io::WRITE_RESULTS );
    src_broker.add_field( src_output_index, source_field );
    src_broker.begin_output_step( src_output_index, 0.0 );
    src_broker.write_defined_output_fields( src_output_index );
    src_broker.end_output_step( src_output_index );


    // TARGET MESH WRITE
    // -----------------

    std::size_t tgt_output_index = tgt_broker.create_output_mesh(
        target_mesh_output_file, stk::io::WRITE_RESULTS );
    tgt_broker.add_field( tgt_output_index, target_field );
    tgt_broker.add_field( tgt_output_index, target_error_field );
    tgt_broker.begin_output_step( tgt_output_index, 0.0 );
    tgt_broker.write_defined_output_fields( tgt_output_index );
    tgt_broker.end_output_step( tgt_output_index );
 


  /*
  // Setup communication.
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
        source_mesh_input_file, "exodus", stk::io::READ_MESH );
    src_broker.set_active_mesh( src_input_index );
    src_broker.create_input_mesh();

    // Add a nodal field to the source part.
    stk::mesh::Field<double>& source_field =
        src_broker.meta_data().declare_field<stk::mesh::Field<double> >(
            stk::topology::NODE_RANK, "u_src" );
    stk::mesh::Part* src_part =
        src_broker.meta_data().get_part( source_mesh_part_name );
    stk::mesh::put_field( source_field, *src_part );

    // Create the source bulk data.
    src_broker.populate_bulk_data();
    Teuchos::RCP<stk::mesh::BulkData> src_bulk_data =
        Teuchos::rcpFromRef( src_broker.bulk_data() );

    // Put some data in the source field. We will use the distance of the node
    // from the origin as the data.
    stk::mesh::Selector src_stk_selector( *src_part );
    stk::mesh::BucketVector src_part_buckets =
        src_stk_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> src_part_nodes;
    stk::mesh::get_selected_entities(
        src_stk_selector, src_part_buckets, src_part_nodes );
    Intrepid::FieldContainer<double> src_node_coords =
        DataTransferKit::STKMeshHelpers::getEntityNodeCoordinates(
            Teuchos::Array<stk::mesh::Entity>(src_part_nodes), *src_bulk_data );
    double* src_field_data;
    int num_src_part_nodes = src_part_nodes.size();
    for ( int n = 0; n < num_src_part_nodes; ++n )
    {
        src_field_data = stk::mesh::field_data( source_field, src_part_nodes[n] );
        src_field_data[0] = dataFunction( src_node_coords(n,0,0),
                                          src_node_coords(n,0,1),
                                          src_node_coords(n,0,2) );
    }


GlobalMPISession mpiSession(&argc,&argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm_default = 
        Teuchos::DefaultComm<int>::getComm();
  int num_procs = comm_default->getSize();

  
  // check input
  int nInCgFile = 0;
  int nOutCgFile = 0;
  
  std::vector<std::string> cgFileName;
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " nCGNSFileIn inCgFileName1 inCgFileName2 ..." << std::endl;
    return 0;
  }
  std::string::size_type sz;   // alias of size_t
  nInCgFile = std::stoi(argv[1],&sz);
  nOutCgFile = nInCgFile; // setting the number of output files the same as input
  for (int iCg=0; iCg < nInCgFile; iCg++)
     cgFileName.push_back(argv[2+iCg]);

  std::cout << "Reading input files #################################################\n";
  // reading cgns file1
  cgnsAnalyzer* cgObj1 = new cgnsAnalyzer(cgFileName[0]);
  cgObj1->loadGrid();
  // attaching partion id to the mesh
  std::vector<double> slnData(cgObj1->getNElement(),0);
  cgObj1->appendSolutionData("partitionOld", slnData, ELEMENTAL, cgObj1->getNElement(), 1);
  // adding new cgns files
  for (int iCg=1; iCg < nInCgFile; iCg++) {
    cgnsAnalyzer* cgObj2 = new cgnsAnalyzer(cgFileName[iCg]);
    cgObj2->loadGrid();
    // appending data
    std::vector<double> slnData(cgObj2->getNElement(),iCg);
    cgObj2->appendSolutionData("partitionOld", slnData, ELEMENTAL, cgObj2->getNElement(), 1);
    // stitching meshes
    cgObj1->stitchMesh(cgObj2, true);
  }
  if (nInCgFile > 1)
     std::cout << "Meshes stitched successfully!\n";

  // checking quality of stitching
  //cgObj1->checkVertex();
  //cgObj1->checkElmConn(3);

  // experiment with mesh partitioner
  //meshPartitioner* mPart1 = new meshPartitioner(cgObj1);
  //mPart1->partition(4);

  std::cout << "Exporting mesh to MAdLib format #####################################\n";
  // exporting mesh to the MAdLib
  MAd::pGModel model = NULL;
  MAd::GM_create(&model,"");
  MAd::pMesh mesh = M_new(model);
  cgObj1->exportToMAdMesh(mesh);
  cgObj1->classifyMAdMeshOpt(mesh);
  
  // writing the mesh to gmsh and convert to vtk
  M_writeMsh(mesh, "stitched.msh", 2, NULL);
  std::cout << "Converting from gmsh to vtk format.\n";
  GModel* trgGModel;
  trgGModel = new GModel("stitched"); 
  trgGModel->readMSH("stitched.msh");
  trgGModel->writeVTK("stitched.vtk", false, true);

  // write physical quantities to vtk file
  std::cout << "Writing physical quantities to vtk file.\n";
  vtkAnalyzer* trgVTK;
  trgVTK = new vtkAnalyzer((char*)"stitched.vtk");
  trgVTK->read();
  trgVTK->report();
  // figure out what is existing on the stitched grid
  int outNData, outNDim;
  std::vector<std::string> slnNameList;
  std::vector<std::string> appSlnNameList;
  cgObj1->getSolutionDataNames(slnNameList);  
  cgObj1->getAppendedSolutionDataName(appSlnNameList);
  slnNameList.insert(slnNameList.end(),
                     appSlnNameList.begin(), appSlnNameList.end());
  // write all data into vtk file
  for (auto is=slnNameList.begin(); is<slnNameList.end(); is++)
  {
    std::vector<double> physData;
    cgObj1->getSolutionDataStitched(*is, physData, outNData, outNDim);
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL)      
      trgVTK->setPointDataArray((*is).c_str(), 1, physData);
    else
      trgVTK->setCellDataArray((*is).c_str(), 1, physData);
  }
  trgVTK->report();
  trgVTK->write("stitchedPhys.vtu");
  */

  /*
  // trying DTK with siera toolkit (stk)
  stk::io::StkMeshIoBroker mesh_data(MPI_COMM_WORLD);
  std::string type = "exodusii";
  std::string file = "./cube_mesh.exo";
  size_t input_index = mesh_data.add_mesh_database(file, type, stk::io::READ_MESH);
  mesh_data.set_active_mesh(input_index);
  mesh_data.create_input_mesh();
  mesh_data.populate_bulk_data(); // expensive step avoid if you can
  */
  
  // query the mesh
  //std::cout << mesh_data.get_coordinate_field().name() << std::endl;
  //std::cout << mesh_data.get_coordinate_field().field_array_rank() << std::endl;
  //mesh_data.get_coordinate_field().get_mesh().dump_all_mesh_info(std::cout);

  // writing to exo format using stk
  // This call adds an output database for results data to mesh_data.
  // No data is written at this time other than verifying that the
  // file can be created on the disk.
  //std::string output_filename = "./out.exo";
  //size_t results_index = mesh_data.create_output_mesh(output_filename, stk::io::WRITE_RESULTS);
  // Iterate all fields and set them as restart fields...
  //const stk::mesh::FieldVector &fields = mesh_data.meta_data().get_fields();
  //for (size_t i=0; i < fields.size(); i++) {
  //  mesh_data.add_field(results_index, *fields[i]); // output
  //}

  // Determine the names of the global fields on the input
  // mesh. These will be used below to define the same fields on the
  // restart and results output databases.
  //std::vector<std::string> global_fields;
  //mesh_data.get_global_variable_names(global_fields);
  //for (auto it=global_fields.begin(); it!=global_fields.end(); it++)
  //  std::cout << *it << std::endl;
  
  /*
  // registering the mesh to DTK
  Teuchos::RCP<DataTransferKit::GeometryManager<DataTransferKit::Box,int> > dtk_srcGeom;

  // Allocate space for 2 boxes
  Teuchos::ArrayRCP<DataTransferKit::Box> boxes( 2 );
  Teuchos::ArrayRCP<int> box_ids( 2 );

  // Create 2 boxes. 1 box will be shared with the other proc.
  boxes[0] = DataTransferKit::Box( 0, 0, 0, 1, 1, 1 );
  box_ids[0] = 0;
  boxes[1] = DataTransferKit::Box( 1, 0, 0, 2, 1, 1 );
  box_ids[1] = 1;
  dtk_srcGeom = Teuchos::rcp(new DataTransferKit::GeometryManager<DataTransferKit::Box,int>(
            boxes, box_ids, comm_default, 3 ));
  */

  /*
  std::cout << "Optimizing the mesh #############################################\n";
  cgObj1->classifyMAdMeshBnd(mesh); // registering boundaries for the optimization step

  // perparing for element quantity interpolation
  std::vector<double> regCntCrdsOld;
  getRegCenters(mesh, regCntCrdsOld);
  //std::vector<double> tmp = getCOG(regCntCrdsOld);
  basicInterpolant* int1 = new basicInterpolant(3, M_numRegions(mesh), 1, regCntCrdsOld);

  // prepare the mesh for the optimization
  std::cout <<"Checking mesh sanity.\n";
  MAd::checkMesh(mesh);
  MAd::M_info(mesh);

  // defining addaptive refinement parameters
  MAd::PWLSField * sizeField = new MAd::PWLSField(mesh);
  sizeField->setCurrentSize();
  sizeField->scale(1.5);
  MAd::MeshAdapter* ma = new MAd::MeshAdapter(mesh,sizeField);
  //ma->uglyTheMesh(0.5,20);
  ma->printParameters();

  // attach all nodal data to the mesh optimizer
  int nData, nDim;
  std::vector<std::string> cgSlnNameList;
  std::vector<std::string> cgAppSlnNameList;
  cgObj1->getSolutionDataNames(cgSlnNameList);  
  cgObj1->getAppendedSolutionDataName(cgAppSlnNameList);
  cgSlnNameList.insert(cgSlnNameList.end(),
                     cgAppSlnNameList.begin(), cgAppSlnNameList.end());
  for (auto is=cgSlnNameList.begin(); is<cgSlnNameList.end(); is++)
  {
    std::vector<double> physData;
    cgObj1->getSolutionDataStitched(*is, physData, nData, nDim);
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL) {
      std::cout << "Reading " << *is << std::endl;
      ma->registerData(*is, physData); 
    }
  }

  // optimize the mesh
  std::cout << "Statistics before optimization: \n";
  ma->printStatistics(std::cout);
  std::cout << "Optimizing the mesh ...\n";
  ma->run();
  std::cout << "Statistics after optimization: \n";
  ma->printStatistics(std::cout);
  // write bulk mesh
  cgObj1->unclassifyMAdMeshBnd(mesh); // remove registered boundaries for proper output
  MAd::M_writeMsh (mesh, "optimizedMesh.msh", 2, NULL);
  std::cout << "Volume Mesh "; 
  MAd::M_info(mesh);
  // convert optimized Gmsh file to vtk
  std::cout << "Converting from optimized gmsh to vtk format.\n";
  trgGModel = new GModel("stitchedOpt"); 
  trgGModel->readMSH("optimizedMesh.msh");
  trgGModel->writeVTK("optimizedMesh.vtk", false, true);

  // skinning the mesh
  std::cout << "Computing surface mesh.\n"; 
  std::vector<int> skinElmIds;
  MAd::pGModel tmpMdl = NULL;
  MAd::GM_create(&tmpMdl,"");
  MAd::pMesh skinMesh = M_new(tmpMdl);
  mesh->skin_me(skinMesh, skinElmIds);
  std::cout << "Num of elements in the mesh = " << MAd::M_numRegions(skinMesh) << std::endl;
  std::cout << "Num of triangles in the mesh = " << MAd::M_numTriangles(skinMesh) << std::endl;
  std::cout << "Minimum surface element index = " << *std::min_element(skinElmIds.begin(),skinElmIds.end()) << std::endl;
  std::cout << "Maximum surface element index = " << *std::max_element(skinElmIds.begin(),skinElmIds.end()) << std::endl;
  skinMesh->classify_unclassified_entities();
  skinMesh->destroyStandAloneEntities();
  //MAd::M_writeMsh (skinMesh, "skinMesh.msh", 2, NULL);
  MAd::SaveGmshMesh (skinMesh, "skinMesh.msh", 2, false, NULL);

  MAd::M_info(skinMesh);

  std::cout << "Transfering solution values #############################################\n";
  // write physical quantities to vtk file
  std::cout << "Writing transfered physical quantities to vtk file.\n";
  // get nodal data after refinement and write them
  vtkAnalyzer* trgVTK2;
  trgVTK2 = new vtkAnalyzer((char*)"optimizedMesh.vtk");
  trgVTK2->read();
  std::cout << "Writing cell and nodal data....\n";
  std::vector<double> regCntCrdsNew;
  getRegCenters(mesh, regCntCrdsNew);
  //std::vector<double> cog = getCOG(regCntCrdsNew);
  int nNewElm = M_numRegions(mesh);
  for (auto is=cgSlnNameList.begin(); is<cgSlnNameList.end(); is++)
  {
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL) {
      std::cout << "Writing nodal " << *is << std::endl;
      std::vector<double> physData;
      ma->getMeshData((*is), &physData);
      trgVTK2->setPointDataArray((*is).c_str(), 1, physData);
      //MAd::NodalDataManagerSgl::instance().writeData((*is),((*is)+".pos").c_str());
    } else {
      //gs field is wiered in irocstar files we dont write it back
      if (!(*is).compare("gs")){
        continue;
      }
      std::cout << "Writing cell-based " << *is << std::endl;
      std::vector<double> oldPhysData;
      std::vector<double> newPhysData;
      int nDataT, nDimT;      
      cgObj1->getSolutionDataStitched(*is, oldPhysData, nDataT, nDimT);
      int1->interpolate(nNewElm, regCntCrdsNew, oldPhysData, newPhysData);
      std::cout << "Size oldPhys = " << oldPhysData.size()
                << " Size newPhys = " << newPhysData.size()
                << " nDataT = " << nDataT
                << std::endl;
      trgVTK2->setCellDataArray((*is).c_str(), 1, newPhysData);
    }
  }
  trgVTK2->report();
  trgVTK2->write("optimizedPhys.vtu");
    
  // partition the mesh
  std::cout << " Partitioning the mesh with METIS.\n";
  //meshPartitioner* mPart = new meshPartitioner(cgObj1);
  meshPartitioner* mPart = new meshPartitioner(mesh);
  mPart->partition(nOutCgFile);
  
  // write partitin ids for the surface mesh
  std::ofstream of;
  of.open("skinMeshPart.dat", std::fstream::trunc);
  std::vector<double> elmPartIds = mPart->getPartedElm();
  for (auto it=skinElmIds.begin(); it!=skinElmIds.end(); it++)
  {
    of << elmPartIds[*it - 1] << "\n";
  } 
  of.close();

  // write CGNS files for the new grid
  for (int iCg=0; iCg<nOutCgFile; iCg++)
  {
     //std::ostringstream exp;
     //exp << "test_" << iCg << ".cgns";
     int1->clearCache();
     std::string fCgName;
     fCgName =cgFileName[iCg];
     std::size_t pos = fCgName.find_last_of("/");
     fCgName = fCgName.substr(pos+1);
     std::cout << "Writing remeshed " << fCgName << std::endl;
     // define elementary information
     cgnsWriter* cgWrtObj = new cgnsWriter(fCgName, cgObj1->getBaseName(), 3, 3);
     cgWrtObj->setUnits(cgObj1->getMassUnit(), cgObj1->getLengthUnit(),
			cgObj1->getTimeUnit(), cgObj1->getTemperatureUnit(),
			cgObj1->getAngleUnit());
     cgWrtObj->setBaseItrData(cgObj1->getBaseItrName(), cgObj1->getNTStep(), cgObj1->getTimeStep());
     cgWrtObj->setZoneItrData(cgObj1->getZoneItrName(), cgObj1->getGridCrdPntr(), cgObj1->getSolutionPntr());
     cgWrtObj->setZone(cgObj1->getZoneName(iCg), cgObj1->getZoneType());
     cgWrtObj->setNVrtx(mPart->getNNdePart(iCg));
     cgWrtObj->setNCell(mPart->getNElmPart(iCg));
     // define coordinates
     cgWrtObj->setGridXYZ(mPart->getCrds(iCg, MAd::M_getVrtXCrds(mesh)), 
			  mPart->getCrds(iCg, MAd::M_getVrtYCrds(mesh)), 
			  mPart->getCrds(iCg, MAd::M_getVrtZCrds(mesh)));
     // define connctivity
     cgWrtObj->setSection(cgObj1->getSectionName(), 
			  (ElementType_t) cgObj1->getElementType(), 
			  mPart->getConns(iCg));
     // define vertex and cell data 
     std::map<std::string, GridLocation_t> slnNLMap = cgObj1->getSolutionNameLocMap();
     for (auto is=slnNLMap.begin(); is!=slnNLMap.end(); is++)
       cgWrtObj->setSolutionNode(is->first, is->second);
     // write skelleton of the file
     cgWrtObj->writeGridToFile();
     // write individual data fields
     std::map<int,std::pair<int,keyValueList> > slnMap = cgObj1->getSolutionMap();
     std::vector<GridLocation_t> gLoc = cgObj1->getSolutionGridLocations();
     std::vector<std::string> slnName = cgObj1->getSolutionNodeNames();
     std::vector<double> regCntCrdsPart = mPart->getElmSlnVec(iCg, regCntCrdsNew, 3);
     //cog = getCOG(regCntCrdsPart);
     int iSol = -1;
     for (auto is=slnMap.begin(); is!=slnMap.end(); is++)
     {
       std::pair<int,keyValueList> slnPair = is->second;
       int slnIdx = slnPair.first;
       keyValueList fldLst = slnPair.second;
       for (auto ifl=fldLst.begin(); ifl!=fldLst.end(); ifl++)
       {
	 iSol++;
	 std::vector<double> stitPhysData;
	 std::vector<double> partPhysData;
	 int nData;
	 if (gLoc[iSol] == Vertex)
	 {
	   nData = mPart->getNNdePart(iCg);
	   ma->getMeshData(ifl->second, &stitPhysData);
	   partPhysData = mPart->getNdeSlnScalar(iCg, stitPhysData);
	 } else {
	   nData = mPart->getNElmPart(iCg);
	   std::vector<double> oldPhysData;
	   int nDataT, nDimT;
	   cgObj1->getSolutionDataStitched(ifl->second, oldPhysData, nDataT, nDimT);
	   int1->interpolate(mPart->getNElmPart(iCg), regCntCrdsPart, oldPhysData, partPhysData);      
           //std::cout << "Minimum element = " 
           //          << *std::min_element(partPhysData.begin(), partPhysData.end())
           //          << "\n Maximum element = " 
           //          << *std::max_element(partPhysData.begin(), partPhysData.end())
           //          << std::endl;
	 }
	 std::cout << "Writing "
		   << nData 
		   << " to "
		   << ifl->second
		   << " located in "
		   << slnName[iSol]
		   << std::endl;
	 // write to file
	 cgWrtObj->writeSolutionField(ifl->second, slnName[iSol], RealDouble, &partPhysData[0]);
       }
     }
     delete cgWrtObj;
  }
  */
  std::cout << "Application ended successfully!\n";
  return 0;
}



///////////////////////////////////////////////////////////////////////
//                       AUX FUNCTIONS                               // 
///////////////////////////////////////////////////////////////////////

/* get element/cell/region center coordinates */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds)
{
   MAd::RIter ri = M_regionIter(msh);
   int rCnt = 0;
   while (MAd::pRegion pr = RIter_next(ri)) 
   {
     double xc[3];
     MAd::R_center(pr, xc);
     regCntCrds.push_back(xc[0]);
     regCntCrds.push_back(xc[1]);
     regCntCrds.push_back(xc[2]);
     /* 
     std::cout << "Region " 
               << rCnt++ 
               << " center coordinate = "
               << xc[0] << " "
               << xc[1] << " "
               << xc[2] << " "
               << std::endl;
     */
   }
} 

/* returns the cartesian coordinates for the geometric center */
std::vector<double> getCOG(std::vector<double>& regCntCrds)
{
  int nNde = regCntCrds.size()/3;
  double x,y,z;
  x=0.0;
  y=0.0;
  z=0.0;
  for (int iNde=0; iNde<nNde; iNde++)
  {
    x += regCntCrds[iNde*3];
    y += regCntCrds[iNde*3 + 1];
    z += regCntCrds[iNde*3 + 2];
  }
  std::vector<double> cog;
  cog.push_back(x/nNde);
  cog.push_back(y/nNde);
  cog.push_back(z/nNde);
  std::cout << "Region center of geometry (" << cog[0]
            << ", " << cog[1] << " , " << cog[2] << ")\n";
  return(cog);
}

//---------------------------------------------------------------------------//
// Data field function.
//---------------------------------------------------------------------------//
double dataFunction( double x, double y, double z )
{
    return std::abs(x) + std::abs(y) + std::abs(z) + 1.0;
}

