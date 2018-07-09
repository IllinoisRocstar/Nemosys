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
#include <vtkMPICommunicator.h>
#include <vtkPointData.h>
#include <vtkGenericCell.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>


#include <mpi.h>
#include <sstream>
#include <cstddef>
#include <set>
#include <AuxiliaryFunctions.H>

vtkStandardNewMacro(GhostGenerator);

GhostGenerator::GhostGenerator()
  : ghostCellMesh(nullptr)
{} 

void GhostGenerator::setPartitions(std::vector<meshBase*>& _partitions)
{
  this->partitions = _partitions;
}

void GhostGenerator::getGlobalIds(int me)
{
  if (this->myGlobalNodeIds.empty())
  {
    this->myGlobToPartNodeMap = partitions[me]->getGlobToPartNodeMap();
    this->myPartToGlobNodeMap = partitions[me]->getPartToGlobNodeMap();
    this->myGlobalNodeIds = getSortedKeys(this->myGlobToPartNodeMap);
      //partitions[me]->getSortedGlobPartNodeIds();
  }
  if (this->myGlobalCellIds.empty())
  {
    this->myGlobalCellIds = getSortedKeys(partitions[me]->getGlobToPartCellMap());
      //partitions[me]->getSortedGlobPartCellIds();
  }
}

void GhostGenerator::getGlobalGhostIds(int me)
{
  if (this->myGlobalGhostNodeIds.empty())
  {
    vtkSmartPointer<vtkIdTypeArray> myVtkGlobalNodeIds
      = vtkIdTypeArray::FastDownCast(
          this->ghostCellMesh->getDataSet()->GetPointData()->GetGlobalIds());
    int start = partitions[me]->getNumberOfPoints();
    this->myGlobalGhostNodeIds.resize(this->ghostCellMesh->getNumberOfPoints()-start);
    for (int i = 0; i < this->myGlobalGhostNodeIds.size(); ++i)
    {
      this->myGlobalGhostNodeIds[i] = (int) myVtkGlobalNodeIds->GetTuple1(i+start);
    }
    std::sort(this->myGlobalGhostNodeIds.begin(), this->myGlobalGhostNodeIds.end());
  }
  if (this->myGlobalGhostCellIds.empty())
  {
    vtkSmartPointer<vtkIdTypeArray> myVtkGlobalCellIds
      = vtkIdTypeArray::FastDownCast(
          this->ghostCellMesh->getDataSet()->GetCellData()->GetGlobalIds());
    int start = partitions[me]->getNumberOfCells();
    this->myGlobalGhostCellIds.resize(this->ghostCellMesh->getNumberOfCells()-start);
    for (int i = 0; i < this->myGlobalGhostCellIds.size(); ++i)
    {
      this->myGlobalGhostCellIds[i] = (int) myVtkGlobalCellIds->GetTuple1(i+start);
    }
    std::sort(this->myGlobalGhostCellIds.begin(), this->myGlobalGhostCellIds.end());
  }
}

void GhostGenerator::getPconnInformation(int me, int numProcs)
{
  this->getGlobalIds(me);
  this->getGlobalGhostIds(me);
  // will hold sent cell's indices
  vtkSmartPointer<vtkIdList> cellIdsList = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  // for each proc that is not me
  for (int i = 0; i < numProcs; ++i)
  {
    if (i == me)
    {
      continue;
    }
    // ------ shared nodes between me and proc i
    std::vector<int> procsGlobalNodeIds(getSortedKeys(partitions[i]->getGlobToPartNodeMap()));
      //partitions[i]->getSortedGlobPartNodeIds());
    // compute the intersection and assign to sharedNodes vec
    std::vector<int> tmpVec;
    std::set_intersection(this->myGlobalNodeIds.begin(), this->myGlobalNodeIds.end(),
                          procsGlobalNodeIds.begin(), procsGlobalNodeIds.end(),
                          std::back_inserter(tmpVec));//(this->sharedNodes[i]));
    this->sharedNodes[i].resize(tmpVec.size());
    // ------ fill shared nodes map and find cells and nodes me sends to proc i
    for (int j = 0; j < tmpVec.size(); ++j)
    {
      cellIdsList->Reset();
      // get local idx of shared node
      int localPntId = myGlobToPartNodeMap[tmpVec[j]];
      // add to sharedNodes with proc i map
      this->sharedNodes[i][j] = localPntId;
      // find cells using this shared node
      partitions[me]->getDataSet()->GetPointCells(localPntId, cellIdsList);
      // for each cell using shared node (these will be cells on the boundary)
      for (int k = 0; k < cellIdsList->GetNumberOfIds(); ++k)
      {
        // get the local cell idx
        int localCellId = cellIdsList->GetId(k);
        // add idx to sentCells to proc i map
        this->sentCells[i].insert(localCellId);
        // get the cell for point extraction
        partitions[me]->getDataSet()->GetCell(localCellId, genCell);
        for (int l = 0; l < genCell->GetNumberOfPoints(); ++l)
        {
          // get node idx of cell point
          int pntCellId = genCell->GetPointId(l);
          // if node is not found in shared nodes, it is sent 
          if (std::find(tmpVec.begin(), tmpVec.end(),
                this->myPartToGlobNodeMap[pntCellId]) == tmpVec.end())
          {
            this->sentNodes[i].insert(pntCellId);
					} 
        }
      }
    }
    // ------ nodes and cells me receives from proc i
    tmpVec.clear();
    std::set_intersection(this->myGlobalGhostNodeIds.begin(), this->myGlobalGhostNodeIds.end(),
                          procsGlobalNodeIds.begin(), procsGlobalNodeIds.end(),
                          std::back_inserter(tmpVec));
    this->receivedNodesNum[i] = tmpVec.size(); 
    tmpVec.clear();
    std::vector<int> procsGlobalCellIds(getSortedKeys(partitions[i]->getGlobToPartCellMap()));
      //partitions[i]->getSortedGlobPartCellIds());
    std::set_intersection(this->myGlobalGhostCellIds.begin(), this->myGlobalGhostCellIds.end(),
                          procsGlobalCellIds.begin(), procsGlobalCellIds.end(),
                          std::back_inserter(tmpVec));
  	this->receivedCellsNum[i] = tmpVec.size();
	}
  std::cout << "bs" << std::endl;
}

void GhostGenerator::writeSharedToPconn(const std::string& type)
{
  auto sharedItr = this->sharedNodes.begin();
  while ( sharedItr != this->sharedNodes.end())
  {
    if (sharedItr->second.size() != 0)
    {
      // num blocks
      pConnVec.push_back(1);
      std::stringstream ss;
      ss << sharedItr->first + 1 << type;
      // zone number
      pConnVec.push_back(std::stoi(ss.str()));
      // number of ids
      pConnVec.push_back(sharedItr->second.size());
      // indices 
      for (int i = 0; i < sharedItr->second.size(); ++i)
      {
        pConnVec.push_back(sharedItr->second[i] + 1);
      }
    }
    ++sharedItr;
  }
}

void GhostGenerator::writeSentToPconn(const std::string& type, bool nodeOrCell)
{
  std::map<int, std::unordered_set<int>>::iterator sentItr;
  std::map<int, std::unordered_set<int>>::iterator end;
  if (!nodeOrCell)
  {
    sentItr = this->sentNodes.begin();
    end = this->sentNodes.end();
  }
  else
  {
    sentItr = this->sentCells.begin();
    end = this->sentCells.end();
  }
  while (sentItr != end)
  {
    if (sentItr->second.size() != 0)
    {
      // num blocks
      pConnVec.push_back(1);
      std::stringstream ss;
      ss << sentItr->first + 1 << type;
      // zone number
      pConnVec.push_back(std::stoi(ss.str()));
      // number of ids
      pConnVec.push_back(sentItr->second.size());
      // indices
  	  auto tmpIt = sentItr->second.begin();
  	  while (tmpIt != sentItr->second.end())
  	  {
        pConnVec.push_back(*tmpIt + 1);
  	  	++tmpIt;
  	  }
    }	
  	++sentItr;
  }
}

void GhostGenerator::writeReceivedToPconn(const std::string& type, int me, bool nodeOrCell)
{

  int offset;
  int end;
  vtkSmartPointer<vtkDataArray> partitionIds;
  if (!nodeOrCell)
  {
    offset = partitions[me]->getNumberOfPoints();
    partitionIds = this->ghostCellMesh->getDataSet()->GetPointData()->GetArray("NodePartitionIds");
    end = ghostCellMesh->getNumberOfPoints();
  }
  else
  {
    offset = partitions[me]->getNumberOfCells();
    partitionIds = this->ghostCellMesh->getDataSet()->GetCellData()->GetArray("CellPartitionIds");
    end = ghostCellMesh->getNumberOfCells();
  }

  while ( offset < end )
  {
    std::cout << "offset " << offset << std::endl;
    double tmp[1];
    partitionIds->GetTuple(offset, tmp);
    int ghostProcId = (int) tmp[0];
    // num blocks
    pConnVec.push_back(1);
    std::stringstream ss;
    ss << ghostProcId + 1 << type;
    // zone number
    pConnVec.push_back(std::stoi(ss.str()));
    std::cout << "zone: " << ss.str() << std::endl;
    std::cout << "ghostprocid : " << ghostProcId+1 << std::endl;
    // number of ids
    int numRecieved = (!nodeOrCell ?
                        this->receivedNodesNum[ghostProcId] :
                        this->receivedCellsNum[ghostProcId]);
    std::cout << "num recieved " << numRecieved << std::endl;
    //std::cout << "offset " << offset << std::endl;
    pConnVec.push_back(numRecieved);
    for (int id = offset; id < numRecieved + offset; ++id)
    {
      pConnVec.push_back(id+1);
    }
    offset += numRecieved;
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
  if (partitions[me]->getNumberOfCells() == 0 || partitions[me] == nullptr)
  {
    if (partitions[me])
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
  ghostGenerator->SetInputData(partitions[me]->getDataSet());
  ghostGenerator->UseGlobalPointIdsOn();
  ghostGenerator->SetGlobalPointIdsArrayName("GlobalNodeIds");
  ghostGenerator->SetHasGlobalCellIds(true);
  ghostGenerator->SetGlobalCellIdsArrayName("GlobalCellIds");
  ghostGenerator->Update();
  std::stringstream ss1;
  ss1 << "ghostProc" << me << ".vtu";
  this->ghostCellMesh 
    = meshBase::Create(vtkDataSet::SafeDownCast(ghostGenerator->GetOutput()),
                                                ss1.str());
  this->ghostCellMesh->write(); 
  this->getPconnInformation(me, numProcs); 
  if (me == 1)
  {

		/*auto it1 = this->receivedNodesNum.begin();
		while (it1 != this->receivedNodesNum.end())
		{
			std::cout << "Num nodes received from proc " << it1->first << " : ";	
			std::cout << it1->second << std::endl;;
			++it1;
		} */	
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
    std::string part = cgObj->getZoneName().substr(0,2);
    std::string type = cgObj->getZoneName().substr(2,3);
    std::cout << "Part and type: " << part << " " << type << std::endl;
    std::cout << "zone as int " << std::stoi(cgObj->getZoneName()) << std::endl;
     
 
    this->writeSharedToPconn(type);
    int notGhostDescriptor = pConnVec.size();
    this->writeSentToPconn(type,0);    
    this->writeReceivedToPconn(type,me,0);
    this->writeSentToPconn(type,1);
    this->writeReceivedToPconn(type,me,1);
    int ghostDescriptor = pConnVec.size() - notGhostDescriptor;
    cgWrtObj->setPconnGhostDescriptor(ghostDescriptor);
    printVec(pConnVec);

    cgWrtObj->setZone(cgObj->getZoneName(), cgObj->getZoneType());
    cgWrtObj->setNVrtx(partitions[me]->getNumberOfPoints());
    cgWrtObj->setNCell(partitions[me]->getNumberOfCells());
    // define coordinates
    std::vector<std::vector<double>> comp_crds(ghostCellMesh->getVertCrds());
    cgWrtObj->setGridXYZ(comp_crds[0], comp_crds[1], comp_crds[2]);
    cgWrtObj->setCoordRind(ghostCellMesh->getNumberOfPoints() - partitions[me]->getNumberOfPoints());
    // define connectivity
    std::vector<int> cgConnReal(partitions[me]->getConnectivities());
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
    cgWrtObj->setNCell(ghostCellMesh->getNumberOfCells() - partitions[me]->getNumberOfCells());
    cgWrtObj->setSection(":T4:virtual", 
                         (ElementType_t) cgObj->getElementType(),
                         cgConnVirtual);
    cgWrtObj->setVirtElmRind(ghostCellMesh->getNumberOfCells() - partitions[me]->getNumberOfCells());
    cgWrtObj->setPconnVec(pConnVec);
    cgWrtObj->writeGridToFile();
    delete cgObj;
    delete cgWrtObj;
  }
}


RocPrepDriver::RocPrepDriver(std::string& fname, int numPartitions)
{
  // ------------- Testing selection of cells by patch no ----------------------------
  
  // loading stitched surf mesh with patch info
  meshBase* stitchedSurf = meshBase::Create("stitchedSurf.vtu");
  // creating selectionId container
  vtkSmartPointer<vtkIdTypeArray> selectionIds = vtkSmartPointer<vtkIdTypeArray>::New();
  selectionIds->SetNumberOfComponents(1);
  // getting patch number array from stitched surf
  vtkSmartPointer<vtkDataArray> patchNumbers 
    = stitchedSurf->getDataSet()->GetCellData()->GetArray("patchNo");
  // looping to pull cell ids with patchNo=5
  for (int i = 0; i < patchNumbers->GetNumberOfTuples(); ++i)
  {
    double tmp[1];
    patchNumbers->GetTuple(i,tmp);
    int patchNo = (int) tmp[0];
    if (patchNo <= 7 && patchNo >= 3)
    {
      selectionIds->InsertNextValue(i);
    }
  }
  // creating a selectionNode which is used to create the selection 
  vtkSmartPointer<vtkSelectionNode> selectionNode
    = vtkSmartPointer<vtkSelectionNode>::New(); 
  selectionNode->SetFieldType(vtkSelectionNode::CELL);
  selectionNode->SetContentType(vtkSelectionNode::INDICES);
  selectionNode->SetSelectionList(selectionIds);
  // create the selection
  vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
  selection->AddNode(selectionNode);
  // creating extractor
  vtkSmartPointer<vtkExtractSelection> extractSelection
    = vtkSmartPointer<vtkExtractSelection>::New();
  // set stitch surf as input on first port
  extractSelection->SetInputData(0, stitchedSurf->getDataSet());
  // set selectionNode as input on second port
  extractSelection->SetInputData(1, selection);
  extractSelection->Update();
  
  meshBase* extractedSurf 
    = meshBase::Create(vtkDataSet::SafeDownCast(
                        extractSelection->GetOutput()), "extracted.vtu");
  extractedSurf->write();
  delete stitchedSurf;
  delete extractedSurf;

  // --------------------------- end testing --------------------------------------


  // load full volume mesh and create METIS partitions
  this->mesh = meshBase::Create(fname);
  this->partitions = meshBase::partition(mesh, numPartitions); 
  // create single process
  vtkSmartPointer<GhostGenerator> p = vtkSmartPointer<GhostGenerator>::New();
  p->setPartitions(this->partitions);
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

RocPrepDriver::~RocPrepDriver()
{
  if (mesh)
  {
    delete mesh;
    mesh = nullptr;
  } 
  if (!partitions.empty())
  {
    for (int i = 0; i < partitions.size(); ++i)
    {
      delete partitions[i];
      partitions[i] = nullptr;
    }
  }
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
}

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t)
{
  std::stringstream names;
  names << case_dir << "/" << prefix << "*" << base_t << "*.cgns";
  return nemAux::glob(names.str());
}
    //auto it = this->sharedNodes.begin();
    //std::cout << "Checking shared nodes...." << std::endl;
    //while ( it != this->sharedNodes.end())
    //{
    //  std::cout << "Nodes shared with proc " << it->first << " : ";
    //  for (int i = 0; i < it->second.size(); ++i)
    //  {
    //    std::cout << it->second[i] << " ";
    //  }
    //  std::cout << std::endl;
    //  ++it;
    //}


		//std::cout << "Checking sent nodes...." << std::endl;
		//auto it2 = this->sentNodes.begin();
		//while (it2 != this->sentNodes.end())
		//{
		//	std::cout << "Nodes sent to proc " << it2->first << " : ";
		//	auto it3 = it2->second.begin();
		//	while (it3 != it2->second.end())
		//	{
		//		std::cout << *it3 <<  " ";
		//		++it3;
		//	}
		//	std::cout << std::endl;	
		//	++it2;
		//}

		//auto it1 = this->receivedNodesNum.begin();
		//while (it1 != this->receivedNodesNum.end())
		//{
		//	std::cout << "Num nodes received from proc " << it1->first << " : ";	
		//	std::cout << it1->second << std::endl;;
		//	++it1;
		//} 	
   
		//std::cout << "Checking sent cells...." << std::endl;
		//it2 = this->sentCells.begin();
		//while (it2 != this->sentCells.end())
		//{
		//	std::cout << "Cells sent to proc " << it2->first << " : ";
		//	auto it3 = it2->second.begin();
		//	while (it3 != it2->second.end())
		//	{
		//		std::cout << *it3 <<  " ";
		//		++it3;
		//	}
		//	std::cout << std::endl;	
		//	++it2;
		//}
		//
		//it1 = this->receivedCellsNum.begin();
		//while (it1 != this->receivedCellsNum.end())
		//{
		//	std::cout << "Num cells received from proc " << it1->first << " : ";
		//	std::cout << it1->second << std::endl;
		//	++it1;
		//}
 
