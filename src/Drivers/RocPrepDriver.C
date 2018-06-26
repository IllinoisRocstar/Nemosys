#include <RocPrepDriver.H>
#include <meshBase.H>
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

#include <mpi.h>
#include <sstream>
#include <cstddef>

// vtkProcess is an abstract class representing a process that can be launched
// by a vtkMultiProcessController. Concrete classes just have to implement
// Execute() method and make sure it set the proper value in ReturnValue.
// The controller has a function consuming vtkProcess objects (SetSingleProcessObject())
//  where the method to execute to is the vtkProcess' Execute() and the data can be 
// encapsulated in the vtkProcess

class MyProcess : public vtkProcess
{
  public:
    static MyProcess* New();
    virtual void Execute();
    void SetPartitions(std::vector<meshBase*>& partitions);
    ~MyProcess() 
    {}

  protected:
    MyProcess();
    std::string volMeshName;
    std::vector<meshBase*> partitions;
    meshBase* dataSet;
    //meshBase* mesh;
};

vtkStandardNewMacro(MyProcess);

MyProcess::MyProcess()
  : dataSet(nullptr)
{} 

void MyProcess::SetPartitions(std::vector<meshBase*>& _partitions)
{
  partitions = _partitions;
}

void MyProcess::Execute()
{
  // setting appropriate return value
  this->ReturnValue = 1;
  // controller isn't set until the controller itself loads MyProcess with
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
  
  // test ghost cell generator
  
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
  //ghostGenerator->Update();
  //ghostGenerator->UpdateInformation();
  //ghostGenerator->UpdatePiece(me,numProcs,1);
  ghostGenerator->Update();

  std::stringstream ss1;
  ss1 << "ghostProc" << me << ".vtu";
  meshBase* ghostCellMesh 
    = meshBase::Create(vtkDataSet::SafeDownCast(ghostGenerator->GetOutput()),
                                                ss1.str());
  ghostCellMesh->write(); 
  delete ghostCellMesh;
}

RocPrepDriver::RocPrepDriver(std::string& fname, int numPartitions)
{
  // load full volume mesh and create METIS partitions
  meshBase* mesh = meshBase::Create(fname);
  std::vector<meshBase*> partitions(meshBase::partition(mesh, numPartitions)); 
  // create single process
  vtkSmartPointer<MyProcess> p = vtkSmartPointer<MyProcess>::New();
  p->SetPartitions(partitions);
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

RocPrepDriver::RocPrepDriver(const std::vector<std::string>& fnames)
{
  std::vector<meshBase*> partitions;
  int numPoints = 0;
  for (int i = 0; i < fnames.size(); ++i)
  {
    partitions.push_back(meshBase::Create(fnames[i]));
    numPoints += partitions[i]->getNumberOfPoints();
  }

  // create global id array
  //vtkSmartPointer<vtkIdTypeArray> globalIds = vtkSmartPointer<vtkIdTypeArray>::New();
  //globalIds->SetName("GlobalNodeIds");
  //globalIds->SetNumberOfValues(numPoints);
  //for (int i = 0; i < numPoints; ++i)
  //{
  //  globalIds->SetValue(i,i);
  //}
   
 
  // create single process
  vtkSmartPointer<MyProcess> p = vtkSmartPointer<MyProcess>::New();
  p->SetPartitions(partitions);
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
  // recreate testing
//  std::vector<std::string> fnames;
//  fnames.push_back(inputjson["Grid 1"].as<std::string>());
//  fnames.push_back(inputjson["Grid 2"].as<std::string>());
//  fnames.push_back(inputjson["Grid 3"].as<std::string>());
//  fnames.push_back(inputjson["Grid 4"].as<std::string>());
//  return new RocPrepDriver(fnames);
}
