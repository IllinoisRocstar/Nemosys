#include <RocPrepDriver.H>
#include <meshBase.H>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkMPIController.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDistributedDataFilter.h>
#include <vtkPUnstructuredGridGhostCellsGenerator.h>
#include <vtkPConnectivityFilter.h>
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
    void SetArgs(int anArgc, char* anArgv[]);
    ~MyProcess() 
    { 
      if (dataSet)
      {
        delete dataSet;
        dataSet = nullptr;
      }
    }

  protected:
    MyProcess();
    int Argc;
    char** Argv;
    meshBase* dataSet;
};

vtkStandardNewMacro(MyProcess);

MyProcess::MyProcess()
  : Argc(0), Argv(nullptr)
{} 

void MyProcess::SetArgs(int anArgc, char* anArgv[])
{
  Argc = anArgc;
  Argv = anArgv; 
}

void MyProcess::Execute()
{
  // setting appropriate return value
  this->ReturnValue = 1;
  // controller isn't set until the controller itself loads MyProcess with
  // SetSingleProcessObject() 
  int numProcs = this->Controller->GetNumberOfProcesses(); 
  int me = this->Controller->GetLocalProcessId();
  // initialize reader for use on procId=0

  /* TODO: removed for testing
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader
    = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  vtkSmartPointer<vtkDataSet> dataSet = nullptr;
  vtkSmartPointer<vtkUnstructuredGrid> unsGrd 
    = vtkSmartPointer<vtkUnstructuredGrid>::New(); */
 
  dataSet = nullptr;
  vtkSmartPointer<vtkUnstructuredGrid> unsGrd 
    = vtkSmartPointer<vtkUnstructuredGrid>::New(); 
  

  int go;
  
  /* TODO: removed for testing
  if (me == 0)
  {
  
    reader->SetFileName(Argv[0]);
    reader->Update();
    dataSet = reader->GetOutput();
    go = 1;
    if ((dataSet->GetNumberOfCells() == 0) || dataSet == nullptr)
    {
      if (dataSet)
      {
        std::cout << "Failure: input file has no cells" << std::endl;
      }
      go = 0;
    }
  } 
  else
  {
    dataSet = unsGrd;
  }*/
  if (me == 0)
  {
    std::string fname(Argv[0]);
    dataSet = meshBase::Create(fname);
    go = 1;
    if ((dataSet->getNumberOfCells() == 0) || dataSet == nullptr)
    {
      if (dataSet)
      {
        std::cerr << "Failure: input file has no cells" << std::endl; 
      }
      go = 0;
    }
  } 
  else
  {
    dataSet = meshBase::Create(unsGrd, "tmp.vtu");
  }

  vtkSmartPointer<vtkMPICommunicator> comm
    = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  comm->Broadcast(&go, 1, 0);
  if (!go)
  {
    return;
  }
  vtkSmartPointer<vtkDistributedDataFilter> ddFilter
    = vtkSmartPointer<vtkDistributedDataFilter>::New();
  // TODO: removed for testing ddFilter->SetInputData(dataSet);
  ddFilter->SetInputData(dataSet->getDataSet());
  ddFilter->SetController(this->Controller);
  ddFilter->UseMinimalMemoryOff();
  //ddFilter->SetBoundaryModeToAssignToAllIntersectingRegions();
  ddFilter->Update();
  ddFilter->SetBoundaryModeToSplitBoundaryCells(); //CLIPS cells at partition interface
  
  // test ghost cell generator
  
  // initialize ghost cell generator
  vtkSmartPointer<vtkPUnstructuredGridGhostCellsGenerator> ghostGenerator
    = vtkSmartPointer<vtkPUnstructuredGridGhostCellsGenerator>::New();
  ghostGenerator->SetController(this->Controller);
  ghostGenerator->SetBuildIfRequired(false);
  ghostGenerator->SetMinimumNumberOfGhostLevels(1);
  ghostGenerator->SetInputConnection(ddFilter->GetOutputPort());
  ghostGenerator->Update();
   
  // get connectivities
  vtkSmartPointer<vtkPConnectivityFilter> connectivity
    = vtkSmartPointer<vtkPConnectivityFilter>::New();
  connectivity->SetInputConnection(ghostGenerator->GetOutputPort());
  connectivity->Update();
 
  // pull out relevant arrays 
  //vtkSmartPointer<vtkDataArray> 

 
  /* TODO: removed for testing
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer1
    = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  std::stringstream ss1;
  ss1 << "ghostProc" << me << ".vtu";
  writer1->SetFileName(&(ss1.str())[0u]);
  writer1->SetInputData(connectivity->GetUnstructuredGridOutput());
  writer1->Write();
 */
  std::stringstream ss1;
  ss1 << "ghostProc" << me << ".vtu";
  meshBase* ghostCellMesh 
    = meshBase::Create(vtkDataSet::SafeDownCast(connectivity->GetUnstructuredGridOutput()),
                                                ss1.str());
  ghostCellMesh->write(); 
  delete ghostCellMesh;
  /* TODO: removed for testing
  // add array for process number for vis
  vtkSmartPointer<vtkUnstructuredGrid> tmp 
    = vtkUnstructuredGrid::SafeDownCast(ddFilter->GetOutput());
  vtkSmartPointer<vtkIntArray> partition = vtkSmartPointer<vtkIntArray>::New();
  partition->SetName("partNo");
  partition->SetNumberOfTuples(tmp->GetNumberOfCells());
  partition->FillComponent(0, me);
  tmp->GetCellData()->AddArray(partition);
  std::cout << "NuM CELLS: " << tmp->GetNumberOfCells() << std::endl;
  // write each partition
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer
    = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  std::stringstream ss;
  ss << "proc" << me << ".vtu";
  std::cout << ss.str() << std::endl; 
  writer->SetFileName(&(ss.str())[0u]);
  writer->SetInputData(tmp);
  writer->Write();
  */
}

RocPrepDriver::RocPrepDriver(std::string& fname)
{
  int retVal = 1;
  vtkSmartPointer<vtkMPIController> controller
    = vtkSmartPointer<vtkMPIController>::New();
  std::vector<char*> arguments;
  arguments.push_back(&(fname[0u]));  
  char** myargv = &arguments[0u];
  int myargc = 2;
  std::cout << *myargv << std::endl;
  controller->Initialize(&myargc,&myargv);
  vtkMultiProcessController::SetGlobalController(controller);
  int me = controller->GetLocalProcessId();
  
  vtkSmartPointer<MyProcess> p = vtkSmartPointer<MyProcess>::New();
  p->SetArgs(myargc,myargv);
  controller->SetSingleProcessObject(p);
  controller->SingleMethodExecute();
  retVal = p->GetReturnValue();
  controller->Finalize();
}

RocPrepDriver* RocPrepDriver::readJSON(json inputjson)
{
  std::string fname = inputjson["Grid To Partition"].as<std::string>();
  return new RocPrepDriver(fname); 
}
