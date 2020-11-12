#include "RefineDriver.H"

#include <fstream>
#include <iostream>

#include "AuxiliaryFunctions.H"
#include "omegahRefineDriver.H"
#include "vtkMesh.H"
#ifdef HAVE_CFMSH
#  include "AMRFoam.H"
#  include "MeshManipulationFoam.H"
#  include "foamMesh.H"
#  include "interpolatePointToCell.H"
#endif
#ifdef MLAMR
#  include <fdeep/fdeep.hpp>
#endif

namespace NEM {
namespace DRV {

// -------------------------------- Refine Driver
// -------------------------------------//
RefineDriver::RefineDriver(const std::string &_mesh, const std::string &method,
                           const std::string &arrayName, double dev_mult,
                           bool maxIsmin, double edgescale,
                           const std::string &ofname, bool transferData,
                           double sizeFactor) {
  std::cout << "RefineDriver created" << std::endl;
  std::cout << "Size Factor = " << sizeFactor << std::endl;
  std::shared_ptr<meshBase> mesh = meshBase::CreateShared(_mesh);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(method, arrayName, dev_mult, maxIsmin, edgescale, ofname,
                   transferData, sizeFactor);
}

RefineDriver::RefineDriver(const std::string &_mesh, const std::string &method,
                           double edgescale, const std::string &ofname,
                           bool transferData) {
  std::cout << "RefineDriver created" << std::endl;
  std::shared_ptr<meshBase> mesh = meshBase::CreateShared(_mesh);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(method, edgescale, ofname, transferData);
}

RefineDriver::RefineDriver(const std::string &_mesh, const std::string &method,
                           const std::string &arrayName, int order,
                           const std::string &ofname, bool transferData) {
  std::shared_ptr<meshBase> mesh = meshBase::CreateShared(_mesh);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(method, arrayName, order, ofname, transferData);
  std::cout << "RefineDriver created" << std::endl;
}

#ifdef HAVE_CFMSH
RefineDriver::RefineDriver(
    const std::string &_mesh, const std::string &ofname,
    const std::string &method, const std::string &inputFile,
    const int &refineInterval, const int &maxRefinement,
    const double &lowerRefineLevel, const double &upperRefineLevel,
    const double &unrefineAbove, const double &unrefineBelow,
    const bool &writeFieldData, const bool &writeMesh,
    const bool &writeRefHistory, const double &timeStep, const double &endTime,
    const int &nBufferLayers, int &maxCells, const std::string &refCriteria,
    const double &startT, const std::string &MLName) {
  std::cout << "RefineDriver created" << std::endl;
  std::shared_ptr<meshBase> mesh = meshBase::CreateShared(_mesh);

  // Initializes AMR workflow
  auto *mshManip = new MeshManipulationFoam();
  mshManip->initAMRWorkFlow(startT, timeStep, endTime);
  if (mshManip) delete mshManip;

  // Starting
  int argc = 1;
  char **argv = new char *[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info << "Create time\n" << Foam::endl;
  Foam::Time runTime(Foam::Time::controlDictName, args);

  int nSteps = endTime / timeStep;
  std::string refCriName = refCriteria;
  nemAux::toLower(refCriName);

  // Reading mesh from vtk file.
  FOAM::foamMesh *fm = new FOAM::foamMesh(mesh);
  fm->write(ofname);

  Foam::polyMesh mesh1(Foam::IOobject(Foam::polyMesh::defaultRegion,
                                      runTime.timeName(), runTime,
                                      Foam::IOobject::MUST_READ));

  // Initialized AMR class
  auto *amr = new Foam::AMRFoam(mesh1);

  if (writeFieldData)
    if (!writeMesh) amr->enableUpdatedField();

  if (writeRefHistory)
    if (!writeMesh) amr->enableRefHistoryData();

  if (writeMesh) amr->enableMeshWriting();

  // Creates scalar field from incoming text/CSV files
  volScalarField meshFieldXY =
      amr->readIncomingCellField(inputFile, "meshFieldXY");

  meshFieldXY.write();

  if (refCriName == "gradient") {
    volScalarField newField = amr->readInitialField("meshFieldXY");
    volScalarField magGrad = amr->getGradient(newField);

    for (int i = 0; i < nSteps; i++) {
      amr->updateAMR(refineInterval, maxRefinement, magGrad, lowerRefineLevel,
                     upperRefineLevel, unrefineAbove, unrefineBelow,
                     nBufferLayers, maxCells);

      runTime++;
    }
    magGrad.write();
  } else if (refCriName == "value") {
    for (int i = 0; i < nSteps; i++) {
      amr->updateAMR(refineInterval, maxRefinement, meshFieldXY,
                     lowerRefineLevel, upperRefineLevel, unrefineAbove,
                     unrefineBelow, nBufferLayers, maxCells);

      runTime++;
    }
    meshFieldXY.write();
#  ifdef MLAMR
  } else if (refCriName == "ml") {
    // Getting ML model input information
    std::vector<double> nonDimUGrad;
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    mesh->getCellDataArray("nonDimUGrad", nonDimUGrad);
    mesh->getCellDataArray("X", X);
    mesh->getCellDataArray("Y", Y);
    mesh->getCellDataArray("Z", Z);

    // Loading ML Model
    const auto model = fdeep::load_model(MLName);

    // Making predictions
    std::vector<int> refinementVec;
    for (int i = 0; i < nonDimUGrad.size(); i++) {
      refinementVec.push_back(0);
    }
    for (int i = 0; i < nonDimUGrad.size(); i++) {
      const auto result = model.predict(
          {fdeep::tensor(fdeep::tensor_shape(static_cast<double>(4)),
                         {X[i], Y[i], Z[i], nonDimUGrad[i]})});

      if (result[0].get(0, 0, 0, 0, 0) >= 0.5) refinementVec[i] = 1;
    }
    volScalarField meshField = amr->assignToVolScalarField(refinementVec);

    for (int i = 0; i < nSteps; i++) {
      amr->updateAMRML(refineInterval, maxRefinement, nBufferLayers, maxCells,
                       meshField);

      runTime++;
    }
#  endif
  } else {
    std::cerr << "Please define refinement operator choice between "
              << "(\"Value\" or \"Gradient\" or \"ML\")"
              << " using \"Refinement Based On\" keyword" << std::endl;
    std::cerr << "ML based refinement can be enabled using ENABLE_MLAMR "
              << "flag!" << std::endl;
    throw;
  }

  amr->writeMesh();

  // Export final mesh to vtu format
  FOAM::foamMesh *fm2 = new FOAM::foamMesh(false);
  fm2->readAMR(runTime);
  vtkMesh *vm = new vtkMesh(fm2->getDataSet(), ofname);
  vm->report();
  vm->write();
  delete vm;

  if (amr) delete amr;

  std::cout << "End!" << std::endl;
}
#endif

RefineDriver::~RefineDriver() {
  std::cout << "RefineDriver destroyed" << std::endl;
}

RefineDriver *RefineDriver::readJSON(const jsoncons::json &inputjson) {
  RefineDriver *refdrvobj;
  std::string _mesh;
  std::string ofname;
  std::string method;
  std::string arrayName;
  double dev_mult;
  bool maxIsmin, transferData;
  double edgescale = 0;
  std::string fieldName = "meshField";
  _mesh = inputjson["Mesh File Options"]["Input Mesh File"].as<std::string>();
  ofname = inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  method =
      inputjson["Refinement Options"]["Refinement Method"].as<std::string>();
  transferData = inputjson["Refinement Options"]["Transfer Data"].as<bool>();
  if (method == "uniform") {
    edgescale = inputjson["Refinement Options"]["Edge Scaling"].as<double>();
    refdrvobj =
        new RefineDriver(_mesh, method, edgescale, ofname, transferData);
  } else if (method == "Z2 Error Estimator") {
    arrayName = inputjson["Refinement Options"]["Array Name"].as<std::string>();
    int order =
        inputjson["Refinement Options"]["Shape Function Order"].as<int>();
    refdrvobj =
        new RefineDriver(_mesh, method, arrayName, order, ofname, transferData);
#ifdef HAVE_CFMSH
  } else if (method == "FV") {
    int refineInterval =
        inputjson["Refinement Options"].contains("Refinement Interval")
            ? inputjson["Refinement Options"]["Refinement Interval"].as<int>()
            : 1;
    int maxRefinement =
        inputjson["Refinement Options"].contains("Maximum Refinement")
            ? inputjson["Refinement Options"]["Maximum Refinement"].as<int>()
            : 1;
    double lowerRefineLevel =
        inputjson["Refinement Options"].contains("Lower Refinement Level")
            ? inputjson["Refinement Options"]["Lower Refinement Level"]
                  .as<double>()
            : -1.0;
    double upperRefineLevel =
        inputjson["Refinement Options"].contains("Upper Refinement Level")
            ? inputjson["Refinement Options"]["Upper Refinement Level"]
                  .as<double>()
            : -1.0;
    double unrefineAbove =
        inputjson["Refinement Options"].contains("Unrefinement Above")
            ? inputjson["Refinement Options"]["Unrefinement Above"].as<double>()
            : -1.0;
    double unrefineBelow =
        inputjson["Refinement Options"].contains("Unrefinement Below")
            ? inputjson["Refinement Options"]["Unrefinement Below"].as<double>()
            : -1.0;
    bool writeFieldData =
        inputjson["Refinement Options"].contains("Write Field Data?")
            ? inputjson["Refinement Options"]["Write Field Data?"].as<bool>()
            : false;
    bool writeMesh =
        inputjson["Refinement Options"].contains("Write Mesh?")
            ? inputjson["Refinement Options"]["Write Mesh?"].as<bool>()
            : false;
    bool writeRefHistory =
        inputjson["Refinement Options"].contains("Write Refinement Data?")
            ? inputjson["Refinement Options"]["Write Refinement Data?"]
                  .as<bool>()
            : false;
    double timeStep =
        inputjson["Refinement Options"].contains("Time Step")
            ? inputjson["Refinement Options"]["Time Step"].as<double>()
            : 1.0;
    double endTime =
        inputjson["Refinement Options"].contains("End Time")
            ? inputjson["Refinement Options"]["End Time"].as<double>()
            : 1.0;
    int nBufferLayers =
        inputjson["Refinement Options"].contains("Buffer Layers")
            ? inputjson["Refinement Options"]["Buffer Layers"].as<int>()
            : 1;
    int maxCells = inputjson["Refinement Options"].contains("Max Cells")
                       ? inputjson["Refinement Options"]["Max Cells"].as<int>()
                       : 500000;
    std::string inFile =
        inputjson["Refinement Options"].contains("Input Field File")
            ? inputjson["Refinement Options"]["Input Field File"]
                  .as<std::string>()
            : "GetError";
    std::string mlmodelname =
        inputjson["Refinement Options"].contains("ML Kernal Name")
            ? inputjson["Refinement Options"]["ML Kernal Name"]
                  .as<std::string>()
            : "fdeep_model.json";
    double startT =
        inputjson["Refinement Options"].contains("Start Time")
            ? inputjson["Refinement Options"]["Start Time"].as<double>()
            : 0.0;
    std::string refCriteria =
        inputjson["Refinement Options"].contains("Refinement Operator")
            ? inputjson["Refinement Options"]["Refinement Operator"]
                  .as<std::string>()
            : "GetError2";

    if (inFile == "GetError") {
      std::cerr << "Please define name of input file using \"Input Field File\""
                << " keyword" << std::endl;
      throw;
    }

    if (refCriteria == "GetError2") {
      std::cerr << "Please define refinement operator "
                << "(\"Value\" or \"Gradient\")"
                << " using \"Refinement Based On\" keyword" << std::endl;
      throw;
    }

    if (lowerRefineLevel == -1 || upperRefineLevel == -1) {
      std::cerr << "Please define lower and upper refinement criteria based On "
                << refCriteria << " of input field using keywords "
                << "\"Lower Refinement Level\" and \"Upper Refinement Level\""
                << std::endl;
      throw;
    }

    refdrvobj = new RefineDriver(
        _mesh, ofname, method, inFile, refineInterval, maxRefinement,
        lowerRefineLevel, upperRefineLevel, unrefineAbove, unrefineBelow,
        writeFieldData, writeMesh, writeRefHistory, timeStep, endTime,
        nBufferLayers, maxCells, refCriteria, startT, mlmodelname);

    return refdrvobj;
#endif
  } else if (method == "Omega_h") {
    refdrvobj = omegahRefineDriver::readJSON(inputjson);
  } else {
    arrayName = inputjson["Refinement Options"]["Array Name"].as<std::string>();
    dev_mult =
        inputjson["Refinement Options"]["StdDev Multiplier"].as<double>();
    maxIsmin =
        inputjson["Refinement Options"]["Max Is Min for Scaling"].as<bool>();
    double sizeFactor;
    sizeFactor =
        inputjson["Refinement Options"].contains("Size Factor")
            ? inputjson["Refinement Options"]["Size Factor"].as<double>()
            : 1.0;
    refdrvobj = new RefineDriver(_mesh, method, arrayName, dev_mult, maxIsmin,
                                 edgescale, ofname, transferData, sizeFactor);
  }
  return refdrvobj;
}

RefineDriver *RefineDriver::readJSON(const std::string &ifname) {
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array()) {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed\n";
    return RefineDriver::readJSON(inputjson[0]);
  } else {
    return RefineDriver::readJSON(inputjson);
  }
}

}  // namespace DRV
}  // namespace NEM
