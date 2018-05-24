#include <meshStitcher.H>
#include <cgnsAnalyzer.H>
#include <meshBase.H>

meshStitcher::meshStitcher(const std::vector<std::string>& _cgFileNames)
  : cgFileNames(_cgFileNames), stitchedMesh(nullptr)
{
  partitions.resize(cgFileNames.size());
  for (int iCg = 0; iCg < cgFileNames.size(); ++iCg)
  {
    partitions[iCg] = new cgnsAnalyzer(cgFileNames[iCg]);
    partitions[iCg]->loadGrid();
    std::vector<double> slnData(partitions[iCg]->getNElement(),iCg);
    partitions[iCg]
      ->appendSolutionData("partitionOld", slnData, ELEMENTAL, partitions[iCg]->getNElement(),1);
    if (iCg)
      partitions[0]->stitchMesh(partitions[iCg],true); 
  } 
  std::cout << "Meshes stitched successfully! #####################################\n";
  std::cout << "Exporting stitched mesh to VTK format #############################\n";
  stitchedMesh = meshBase::Create(partitions[0]->getVTKMesh(),"stitched.vtu");
  std::cout << "Transferring physical quantities to vtk mesh ######################\n";
  // figure out what is existing on the stitched grid
  int outNData, outNDim;
  std::vector<std::string> slnNameList;
  std::vector<std::string> appSlnNameList;
  partitions[0]->getSolutionDataNames(slnNameList);  
  partitions[0]->getAppendedSolutionDataName(appSlnNameList);
  slnNameList.insert(slnNameList.end(),
                     appSlnNameList.begin(), appSlnNameList.end());
  // write all data into vtk file
  for (auto is=slnNameList.begin(); is<slnNameList.end(); is++)
  {
    std::vector<double> physData;
    partitions[0]->getSolutionDataStitched(*is, physData, outNData, outNDim);
    solution_type_t dt = partitions[0]->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL)  
    {      
      std::cout << "Writing nodal " << *is << std::endl; 
      stitchedMesh->setPointDataArray((*is).c_str(), physData);
    }
    else
    {
      // gs field is 'weird' in irocstar files- we don't write it back
      if (!(*is).compare("gs"))
        continue;
      std::cout << "Writing cell-based " << *is << std::endl;
      stitchedMesh->setCellDataArray((*is).c_str(), physData);
    }
  }
  stitchedMesh->report();
  stitchedMesh->write();
}

meshStitcher::~meshStitcher()
{
  for (int i = 0; i < partitions.size(); ++i)
  {    
    if (partitions[i])
    {
      delete partitions[i];
      partitions[i] = 0;
    }
  }
  if (stitchedMesh)
  {
    delete stitchedMesh;
    stitchedMesh = nullptr;
  }
}

meshBase* meshStitcher::getStitchedMB()
{
  if (stitchedMesh)
    return stitchedMesh;
  else
  {
    std::cerr << "No stitched mesh to return!" << std::endl;
    exit(1);
  }
}

cgnsAnalyzer* meshStitcher::getStitchedCGNS()
{
  if (partitions[0])
    return partitions[0];
  else
  {
    std::cerr << "No stitched mesh to return!" << std::endl;
    exit(1);
  }
}
