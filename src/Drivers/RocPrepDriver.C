#include <RocPrepDriver.H>
#include <meshBase.H>
#include <cgnsWriter.H>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkCellData.h>
#include <vtkMPIController.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDistributedDataFilter.h>
#include <vtkPUnstructuredGridGhostCellsGenerator.h>
#include <vtkPConnectivityFilter.h>
#include <vtkRemoveGhosts.h>
#include <vtkMPICommunicator.h>
#include <vtkProcess.h>
#include <vtkPointData.h>

#include <mpi.h>
#include <sstream>
#include <cstddef>
#include <set>

// vtkProcess is an abstract class representing a process that can be launched
// by a vtkMultiProcessController. Concrete classes just have to implement
// Execute() method and make sure it set the proper value in ReturnValue.
// The controller has a function consuming vtkProcess objects (SetSingleProcessObject())
//  where the method to execute to is the vtkProcess' Execute() and the data can be 
// encapsulated in the vtkProcess

class GhostGenerator : public vtkProcess
{
  public:
    static GhostGenerator* New();
    virtual void Execute();
    void setPartitions(std::vector<meshBase*>& partitions);
    void getSharedNodes(int me, int numProcs);

    ~GhostGenerator() 
    {
      delete dataSet;
    }

  protected:
    GhostGenerator();
    std::vector<meshBase*> partitions;
    meshBase* dataSet;
    std::map<int, std::vector<int>> sharedNodes;
    std::vector<int> sendNodes;
    std::vector<int> recievedNodes;
    std::vector<int> sentCells;
    std::vector<int> recievedCells;
 
};

vtkStandardNewMacro(GhostGenerator);

GhostGenerator::GhostGenerator()
  : dataSet(nullptr)
{} 

void GhostGenerator::setPartitions(std::vector<meshBase*>& _partitions)
{
  partitions = _partitions;
}

void GhostGenerator::getSharedNodes(int me, int numProcs)
{
  // get global node array of me
  vtkSmartPointer<vtkIdTypeArray> myVtkGlobalNodeIds
    = vtkIdTypeArray::FastDownCast(
        partitions[me]->getDataSet()->GetPointData()->GetGlobalIds());
  // convert to vector and sort
  std::vector<int> myGlobalNodeIds(myVtkGlobalNodeIds->GetNumberOfTuples());
  for (int i = 0; i < myGlobalNodeIds.size(); ++i)
  {
    myGlobalNodeIds[i] = (int) myVtkGlobalNodeIds->GetTuple1(i);
  }
  std::sort(myGlobalNodeIds.begin(), myGlobalNodeIds.end());

  // for each proc that is not me
  for (int i = 0; i < numProcs; ++i)
  {
    if (i == me)
    {
      continue;
    }
    // get procs global node ids
    vtkSmartPointer<vtkIdTypeArray> procsVtkGlobalNodeIds
      = vtkIdTypeArray::FastDownCast(
          partitions[i]->getDataSet()->GetPointData()->GetGlobalIds());
    // convert to vector and sort
    std::vector<int> procsGlobalNodeIds(procsVtkGlobalNodeIds->GetNumberOfTuples());
    for (int j = 0; j < procsGlobalNodeIds.size(); ++j)
    {
      procsGlobalNodeIds[j] = (int) procsVtkGlobalNodeIds->GetTuple1(j);
    }
    std::sort(procsGlobalNodeIds.begin(), procsGlobalNodeIds.end());
    // compute the intersection and assign to sharedNodes vec
    std::vector<int> tmpShared;
    std::set_intersection(myGlobalNodeIds.begin(), myGlobalNodeIds.end(),
                          procsGlobalNodeIds.begin(), procsGlobalNodeIds.end(),
                          std::back_inserter(this->sharedNodes[i]));
  }
}

void GhostGenerator::Execute()
{
  // setting appropriate return value
  this->ReturnValue = 1;
  // controller isn't set until the controller itself loads GhostGenerator with
  // SetSingleProcessObject() 
  int numProcs = this->Controller->GetNumberOfProcesses(); 
  if (numProcs != partitions.size())
  {
    std::cerr << "mpirun was called with a different number of processors than requested"
              << " partitions.\n These need to be equal!" << std::endl;
    this->ReturnValue = 0;
    return;
  } 
  int me = this->Controller->GetLocalProcessId();
  int go;
  // load partition per process
  dataSet = partitions[me];
  if (dataSet->getNumberOfCells() == 0 || dataSet == nullptr)
  {
    if (dataSet)
    {
      std::cerr << "Failure: input mesh has no cells" << std::endl;
    }
    go = 0;
  }
  go = 1;
  vtkSmartPointer<vtkMPICommunicator> comm
    = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  comm->Broadcast(&go, 1, 0);
  if (!go)
  {
    this->ReturnValue = 0;
    return;
  }
  
  // initialize ghost cell generator
  vtkSmartPointer<vtkPUnstructuredGridGhostCellsGenerator> ghostGenerator
    = vtkSmartPointer<vtkPUnstructuredGridGhostCellsGenerator>::New();
  ghostGenerator->SetController(this->Controller);
  ghostGenerator->BuildIfRequiredOff();
  ghostGenerator->SetMinimumNumberOfGhostLevels(1);
  ghostGenerator->SetInputData(dataSet->getDataSet());
  ghostGenerator->UseGlobalPointIdsOn();
  ghostGenerator->SetGlobalPointIdsArrayName("GlobalNodeIds");
  ghostGenerator->SetHasGlobalCellIds(true);
  ghostGenerator->SetGlobalCellIdsArrayName("GlobalCellIds");
  ghostGenerator->Update();
  std::stringstream ss1;
  ss1 << "ghostProc" << me << ".vtu";
  meshBase* ghostCellMesh 
    = meshBase::Create(vtkDataSet::SafeDownCast(ghostGenerator->GetOutput()),
                                                ss1.str());
  ghostCellMesh->write(); 
  
  this->getSharedNodes(me, numProcs); 
 
  if (me == 0)
  {

    std::cout << "Checking shared nodes...." << std::endl;
    auto it = this->sharedNodes.begin();
    while ( it != this->sharedNodes.end())
    {
      for (int i = 0; i < it->second.size(); ++i)
      {
        std::cout << it->second[i] << " ";
      }
      std::cout << std::endl;
    }
    
    
    // write cgns file
    
    cgnsAnalyzer* cgObj = new cgnsAnalyzer("fluid-grid_00.000000_00001.cgns");
    cgObj->loadGrid(1);
    std::string fCgName("fluid-grid_00.000000_00001New.cgns");
    cgnsWriter* cgWrtObj = new cgnsWriter(fCgName, cgObj->getBaseName(), 3, 3);
    // define elementary information
    cgWrtObj->setUnits(cgObj->getMassUnit(), cgObj->getLengthUnit(),
                       cgObj->getTimeUnit(), cgObj->getTemperatureUnit(),
                       cgObj->getAngleUnit());
    cgWrtObj->setBaseItrData(cgObj->getBaseItrName(),
                             cgObj->getNTStep(),
                             cgObj->getTimeStep());
    cgWrtObj->setZoneItrData(cgObj->getZoneItrName(),
                             cgObj->getGridCrdPntr(),
                             cgObj->getSolutionPntr());
    cgWrtObj->setZone(cgObj->getZoneName(), cgObj->getZoneType());
    cgWrtObj->setNVrtx(dataSet->getNumberOfPoints());
    cgWrtObj->setNCell(dataSet->getNumberOfCells());
    // define coordinates
    std::vector<std::vector<double>> comp_crds(ghostCellMesh->getVertCrds());
    cgWrtObj->setGridXYZ(comp_crds[0], comp_crds[1], comp_crds[2]);
    cgWrtObj->setCoordRind(ghostCellMesh->getNumberOfPoints() - dataSet->getNumberOfPoints());
    // define connectivity
    std::vector<int> cgConnReal(dataSet->getConnectivities());
    for (auto it = cgConnReal.begin(); it != cgConnReal.end(); ++it)
    {
      *it += 1;
    }
    std::vector<int> cgConnVirtual(ghostCellMesh->getConnectivities());
    cgConnVirtual.erase(cgConnVirtual.begin(), cgConnVirtual.begin()+cgConnReal.size());
    for (auto it = cgConnVirtual.begin(); it != cgConnVirtual.end(); ++it)
    {
      *it += 1;
    }
    cgWrtObj->setSection(cgObj->getSectionName(), 
                         (ElementType_t) cgObj->getElementType(),
                         cgConnReal);
    cgWrtObj->setNCell(ghostCellMesh->getNumberOfCells() - dataSet->getNumberOfCells());
    cgWrtObj->setSection(":T4:virtual", 
                         (ElementType_t) cgObj->getElementType(),
                         cgConnVirtual);
    cgWrtObj->writeGridToFile();
    delete cgObj;
    delete cgWrtObj;
  }
  delete ghostCellMesh;
}


RocPrepDriver::RocPrepDriver(std::string& fname, int numPartitions)
{
  // load full volume mesh and create METIS partitions
  meshBase* mesh = meshBase::Create(fname);
  std::vector<meshBase*> partitions(meshBase::partition(mesh, numPartitions)); 
  // create single process
  vtkSmartPointer<GhostGenerator> p = vtkSmartPointer<GhostGenerator>::New();
  p->setPartitions(partitions);
  // intialize mpi
  MPI_Init(NULL,NULL);
  // initialize controller 
  vtkSmartPointer<vtkMPIController> controller
    = vtkSmartPointer<vtkMPIController>::New();
  controller->Initialize();
  vtkMultiProcessController::SetGlobalController(controller);
  // pass single process to controller
  controller->SetSingleProcessObject(p);
  controller->SingleMethodExecute();
  if(!p->GetReturnValue())
  {
    exit(1);
  }
  controller->Finalize();
  delete mesh;
}

RocPrepDriver::RocPrepDriver(const std::vector<std::string>& fnames)
{
  std::vector<meshBase*> partitions;
  int numPoints = 0;
  for (int i = 0; i < fnames.size(); ++i)
  {
    partitions.push_back(meshBase::Create(fnames[i]));
    numPoints += partitions[i]->getNumberOfPoints();
  }

  // create single process
  vtkSmartPointer<GhostGenerator> p = vtkSmartPointer<GhostGenerator>::New();
  p->setPartitions(partitions);
  // intialize mpi
  MPI_Init(NULL,NULL);
  // initialize controller 
  vtkSmartPointer<vtkMPIController> controller
    = vtkSmartPointer<vtkMPIController>::New();
  controller->Initialize();
  vtkMultiProcessController::SetGlobalController(controller);
  // pass single process to controller
  controller->SetSingleProcessObject(p);
  controller->SingleMethodExecute();
  if(!p->GetReturnValue())
  {
    exit(1);
  }
  controller->Finalize();
}

RocPrepDriver* RocPrepDriver::readJSON(json inputjson)
{
  std::string fname = inputjson["Grid To Partition"].as<std::string>();
  int numPartitions = inputjson["Number Of Partitions"].as<int>();
  return new RocPrepDriver(fname, numPartitions);
  //std::string case_dir = inputjson["Case Directory"].as<std::string>();
  //std::string base_t = inputjson["Base Time Step"].as<std::string>();
  //json templates = inputjson["Template CGNS Files"];
  //std::string templateCgnsName = inputjson["Template CGNS"].as<std::string>(); 
  //bool readAll = false;
  //std::vector<std::string> prefixes;
  //int i = 1;
  //while (!readAll)
  //{
  //  std::string prefix_key = "File Prefix";
  //  prefix_key += std::to_string(i);  
  //  if(templateCgnsName.has_key(prefix_key))
  //    prefixes.push_back(templateCgnsName[prefix_key].as<std::string>();  
  //  else
  //    readAll = true;
  //  ++i;
  //}
  //return new RocPrepDriver(case_dir, base_t, prefixes, fname, numPartitions); 
  // recreate testing
//  std::vector<std::string> fnames;
//  fnames.push_back(inputjson["Grid 1"].as<std::string>());
//  fnames.push_back(inputjson["Grid 2"].as<std::string>());
//  fnames.push_back(inputjson["Grid 3"].as<std::string>());
//  fnames.push_back(inputjson["Grid 4"].as<std::string>());
//  return new RocPrepDriver(fnames);
}

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
}
