/* implementation of mesh partition class(es) */

#include <meshPartitioner.H>
#include <cgnsAnalyzer.H>
#include <meshBase.H>

/* Implementation of meshPartition class */
//meshPartition::meshPartition(int pidx, std::vector<int> glbNdePartedIdx, std::vector<int> glbElmPartedIdx)
meshPartition::meshPartition(int pidx, std::vector<int> glbElmPartedIdx, std::vector<int> glbElmConn, MeshType_t inMshType)
{
  mshType = inMshType;
  pIdx = pidx;
  nNde = 0;
  nElm = 0;
  if (mshType == MESH_TETRA_4) 
    nNdeElm = 4;
  else if (mshType == MESH_TRI_3)
    nNdeElm = 3; 
  // finding partition elements
  int glbElmIdx = 1;
  for (auto ie=glbElmPartedIdx.begin(); ie!=glbElmPartedIdx.end(); ie++)
  {
    if (*ie==pIdx)
    {
      nElm++;
      elmIdxGlobToPart[glbElmIdx] = nElm;
      elmIdxPartToGlob[nElm] = glbElmIdx;
      globElmIdx.push_back(glbElmIdx);
      //std::cout << "My element = " << glbElmIdx << std::endl;
      // attaching nodes of this element to the current partition
      for (int iNde=0; iNde<nNdeElm; iNde++)
      {
        int glbNdeIdx = glbElmConn[(glbElmIdx-1)*nNdeElm + iNde];
        auto it = ndeIdxGlobToPart.find(glbNdeIdx);
        int partNdeIdx;
        if (it==ndeIdxGlobToPart.end())
        {
          nNde++;
          partNdeIdx = nNde;
          ndeIdxGlobToPart[glbNdeIdx] = nNde;
          ndeIdxPartToGlob[nNde] = glbNdeIdx;
          globNdeIdx.push_back(glbNdeIdx);
          //std::cout << "partition " << pIdx
          //          << " global node id " << glbNdeIdx
          //          << " is local " << nNde
          //          << std::endl;
        } else {
          partNdeIdx = it->second;
        }
        partElmConn.push_back(partNdeIdx);
      }
    }
    glbElmIdx++;
  }
  //std::cout <<  "Min connectivity passed to partitioner = "
  //          << *std::min_element(glbElmConn.begin(), glbElmConn.end()) << std::endl;

  // Debug information
  //std::cout << "I am partition " << pidx 
  //          << "\n cycled on " << glbElmIdx-1 << " global elements "
  //          << "\n my min global element idx is " << *std::min_element(globElmIdx.begin(), globElmIdx.end())
  //          << "\n my max global element idx is " << *std::max_element(globElmIdx.begin(), globElmIdx.end())
  //          << "\n my nNde = " << nNde
  //          << " my nElm = " << nElm << std::endl;
}

std::vector<double> meshPartition::getCrds(std::vector<double>& crds)
{
  std::vector<double> x;
  for (auto in=globNdeIdx.begin(); in!=globNdeIdx.end(); in++)
    x.push_back(crds[*in-1]);
  return(x);
}

std::vector<int> meshPartition::getConns()
{
  return(partElmConn);
}


std::vector<double> meshPartition::getElmSlns(std::vector<double>& slns)
{
  std::vector<double> x;
  for (auto ie=globElmIdx.begin(); ie!=globElmIdx.end(); ie++)
    x.push_back(slns[*ie-1]);
  return(x);
}

std::vector<double> meshPartition::getElmSlnsVec(std::vector<double>& slns, int nComp)
{
  std::vector<double> x;
  for (auto ie=globElmIdx.begin(); ie!=globElmIdx.end(); ie++)
    for (int iComp=0; iComp<nComp; iComp++)
      x.push_back(slns[(*ie-1)*nComp + iComp]);
  return(x);
}

std::map<int,int> meshPartition::getPartToGlobNodeMap()
{
  return ndeIdxPartToGlob;
}  

std::map<int,int> meshPartition::getPartToGlobElmMap()
{
  return elmIdxPartToGlob;
}  

/* Implementation of class partitioner class */
meshPartitioner::meshPartitioner(cgnsAnalyzer* inCg)
{
  nNde = inCg->getNVertex();
  nElm = inCg->getNElement();
  elmConnVec = inCg->getElementConnectivity(-1);
  elmConn.insert(elmConn.begin(), elmConnVec.begin(), elmConnVec.end());
  // 0-indexing connectivity 
  for (auto it=elmConn.begin(); it!=elmConn.end(); it++)
    *it = *it-1;
  nPart = 0;
  // coverting between CGNS to local type
  switch (inCg->getElementType())
  {
    case CG_TETRA_4:
      meshType = MESH_TETRA_4;
      break;
    case CG_TRI_3:
      meshType = MESH_TRI_3;
      break;
    default:
      std::cerr << "Unknown element type!\n";
      break;
  }
}


meshPartitioner::meshPartitioner(const meshBase* inMB)
{
  vtkSmartPointer<vtkCellTypes> celltypes = vtkSmartPointer<vtkCellTypes>::New();
  inMB->getDataSet()->GetCellTypes(celltypes);
  for (int i = 0; i < celltypes->GetNumberOfTypes(); ++i)
  {
    if (!(celltypes->IsType(VTK_TETRA)) && !(celltypes->IsType(VTK_TRIANGLE)))
    {
      std::cerr << "Partitioner works for 4-node tet and 3-node tri elements only.\n";
      exit(-1);
    }
  }
  nNde = inMB->getNumberOfPoints();
  elmConnVec = inMB->getConnectivities(); 
  elmConn = elmConnVec;
  // 1 based idex for elmConnVec, 0 based for elmConn
  for (int i= 0; i < elmConnVec.size(); ++i)
  {
    elmConnVec[i] = elmConnVec[i] + 1;
  }
  nElm = inMB->getNumberOfCells();
  if (celltypes->IsType(VTK_TETRA))
    meshType = MESH_TETRA_4;
  else if (celltypes->IsType(VTK_TRIANGLE))
    meshType = MESH_TRI_3;
  else
  {
    std::cerr << "Mesh with unspported element types.\n";
    exit(-1);
  }
  nPart = 0;
  std::cout << " ------------------- Partitioner Stats ---------------------\n";
  std::cout << "Mesh type : " << ( meshType==0 ? "Triangular":"Tetrahedral") << std::endl;
  std::cout << "Number of nodes = " << nNde
            << "\nNumber of elements = " << nElm << std::endl;
  std::cout << "Size of elmConn = " << elmConnVec.size() << std::endl;
  std::cout <<  "Min connectivity passed to partitioner = "
            << *std::min_element(elmConnVec.begin(), elmConnVec.end()) << std::endl;
  std::cout <<  "Max connectivity passed to partitioner = "
            << *std::max_element(elmConnVec.begin(), elmConnVec.end()) << std::endl;
  std::cout << " ----------------------------------------------------------\n";
}


meshPartitioner::meshPartitioner(MAd::pMesh inMesh)
{
  //  only implemented for TETRA_4 elements 
  if (inMesh->nbQuads!=0 || inMesh->nbHexes!=0 || inMesh->nbPrisms!=0)
  {
    std::cerr << "Partitioner works 4-node tet and 3-node tri elements only.\n";
    exit(-1);
  }
  // during adaptation a series of new nodes will be created but
  // MAdLib does not destroy old stand alone points so have to use 
  // maxId instead of number of nodes.
  nNde = MAd::M_numVertices(inMesh); 
  elmConnVec = MAd::M_getConnectivities(inMesh);
  elmConn.insert(elmConn.begin(), elmConnVec.begin(), elmConnVec.end());     
  // 0-indexing connectivity 
  for (auto it=elmConn.begin(); it!=elmConn.end(); it++)
    *it = *it-1;
  if (M_numTets(inMesh) > 0)
  {
    nElm = MAd::M_numRegions(inMesh);
    meshType = MESH_TETRA_4;
  } else if (M_numTriangles(inMesh) > 0) {
    nElm = MAd::M_numTriangles(inMesh);
    meshType = MESH_TRI_3;
  } else {
    std::cerr << "Mesh with unsupported element types.\n";
    exit(-1);
  }
  nPart = 0;
  std::cout << " ------------------- Partitioner Stats ---------------------\n";
  std::cout << "Mesh type : " << ( meshType==0 ? "Triangular":"Tetrahedral") << std::endl;
  std::cout << "Number of nodes = " << nNde
            << "\nNumber of elements = " << nElm << std::endl;
  std::cout << "Size of elmConn = " << elmConnVec.size() << std::endl;
  std::cout <<  "Min connectivity passed to partitioner = "
            << *std::min_element(elmConnVec.begin(), elmConnVec.end()) << std::endl;
  std::cout <<  "Max connectivity passed to partitioner = "
            << *std::max_element(elmConnVec.begin(), elmConnVec.end()) << std::endl;
  std::cout << " ----------------------------------------------------------\n";
}

void meshPartitioner::setNPartition(int nPartition)
{
  nPart = nPartition;
}

int meshPartitioner::getNPartition()
{
  return(nPart);
}

int meshPartitioner::partition(int nPartition)
{
  setNPartition(nPartition);
  return(partition() );
}

int meshPartitioner::partition()
{
  // check
  if (nPart == 0)
  {
    std::cerr << "Number of partitions must be set.\n";
    exit(-1);
  }

  std::cout << "Partitioning the mesh.\n";
  // prepare metis datastructs
  idx_t eptr[nElm+1];
  idx_t objval = 0;
  epart.resize(nElm,0);
  npart.resize(nNde,0);
  int ncommon = 1;
  switch (meshType) {
     case MESH_TETRA_4:
       eptr[0] = 0;
       for (int iElm=1; iElm<=nElm; iElm++)
   eptr[iElm] = iElm*4;
       ncommon = 3;
       break;
     case MESH_TRI_3:
       eptr[0] = 0;
       for (int iElm=1; iElm<=nElm; iElm++)
   eptr[iElm] = iElm*3;
       ncommon = 2;
       break;
     default:
       std::cerr << "Unknown or unimplemented element type.\n";    
  }
  // setting options (some default values, should be tailored)
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0; // 0-based index
  options[METIS_OPTION_CONTIG] = 1;    // try for contiguous partitions
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // multilevel k-way partitioning
  //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB; // multileve recursive bisectioning
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // minimize edge cut METIS_OBJTYPE_VOL(comm vol)
  //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_UFACTOR] = 1; // load imbalance of 1.001 
  //options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // sorted heavy-edge matching

  // To try to match METIS in Rocstar partitioner, the options can be changed to:
  // - Add SHEM method:      options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
  // - Replace KWAY with RB: options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB
  // - Replace VOL with CUT: options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT

  //options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME | 
  //                               METIS_DBG_CONNINFO | METIS_DBG_CONTIGINFO;   
  int res;

  /*
  // DEBUG: checking equivalent graph of the mesh
  std::cout << "METIS: Calculating mesh dual graph...." << std::endl;
  idx_t *xadj, *adjncy;
  idx_t numflag = 0;
  res = METIS_MeshToDual(&nElm, &nNde, &eptr[0], &elmConn[0], &ncommon, &numflag, 
                   &xadj, &adjncy);
  std::ofstream graphFile;
  graphFile.open("dualGraph.txt");
  for (int iElm=0; iElm<nElm; iElm++)
  {
    graphFile << " Elm Idx " << iElm+1 << " -> ";
    for (int iAdj=xadj[iElm]; iAdj<xadj[iElm+1]; iAdj++)
      graphFile << adjncy[iAdj] << " ";
    graphFile << " ***  " << std::endl;
  }
  graphFile.close();
  std::cout << "Dual graph for the mesh is written to -> dualGraph.txt " << std::endl;
  */

  // calling metis partioner
  res = METIS_OK;
  if (nPart > 1)
  {
    /*
    // DEBUG Information
    std::cout << "Sending data to METIS..." << std::endl;
    std::cout << "nElm = " << nElm << std::endl;
    std::cout << "nNde = " << nNde << std::endl;
    std::cout << "eptr[0] = " << eptr[0] << std::endl;
    std::cout << "eptr[end] = " << eptr[nElm] << std::endl;
    std::cout << "eind = " << elmConn[0] << std::endl;
    */
  
   // res = METIS_PartMeshNodal(&nElm, &nNde, &eptr[0], &elmConn[0], NULL, NULL,
   //                           &nPart, NULL, options, &objval, &epart[0], &npart[0]);
    res = METIS_PartMeshDual(&nElm, &nNde, &eptr[0], &elmConn[0], NULL, NULL, 
                             &ncommon, &nPart, NULL, options, &objval, &epart[0], &npart[0]);
    
    std::cout << "Received data from METIS" << std::endl;
  }
  // output partitioning results
  // check success
  if (res == METIS_OK)
  {
    std::cout << "Successfully partitioned the mesh.\n";
    // removing the first member of the npart as default 
    // index starts from 1
    npart.erase(npart.begin());
    buildPartitions();
    return 0;
  } 
  else 
  {
    std::cerr << "Failed to partition with error code " << res << std::endl;
    return 1;
  }
}


std::vector<double> meshPartitioner::getPartedNde()
{
  std::vector<double> ndeParted(npart.begin(), npart.end());
  return(ndeParted);
}

std::vector<double> meshPartitioner::getPartedElm()
{
  std::vector<double> elmParted(epart.begin(), epart.end());
  return(elmParted);
}

void meshPartitioner::setPartedElm(std::vector<double>& prtElm)
{
  int minp, maxp;
  minp = *std::min_element(prtElm.begin(), prtElm.end());
  maxp = *std::max_element(prtElm.begin(), prtElm.end());
  setNPartition(maxp-minp+1);
  epart.insert(epart.begin(), prtElm.begin(), prtElm.end());
  buildPartitions();
}

void meshPartitioner::buildPartitions()
{
  for (int iPart=0; iPart<nPart; iPart++)
  {
    meshPartition* newPart = new meshPartition(iPart, epart, elmConnVec, meshType);
    meshParts.push_back(newPart);
  }
}

std::vector<double> meshPartitioner::getCrds(int iPart, std::vector<double> crds)
{
  return(meshParts[iPart]->getCrds(crds));
}

std::vector<int> meshPartitioner::getConns(int iPart)
{
  return(meshParts[iPart]->getConns());
}

int meshPartitioner::getNNdePart(int iPart)
{
  return(meshParts[iPart]->nNde);
}

int meshPartitioner::getNElmPart(int iPart)
{
  return(meshParts[iPart]->nElm);
}

std::vector<double> meshPartitioner::getNdeSlnScalar(int iPart, std::vector<double>& slns)
{
  // for scalar nodal values the same operator acting on 
  // coordinates can be used to find solution values for
  // the partition
  return(getCrds(iPart, slns));
}

std::vector<double> meshPartitioner::getElmSlnScalar(int iPart, std::vector<double>& slns)
{
  return(meshParts[iPart]->getElmSlns(slns));
}

std::vector<double> meshPartitioner::getElmSlnVec(int iPart, std::vector<double>& slns, int nComp)
{
  return(meshParts[iPart]->getElmSlnsVec(slns, nComp));
}
  

std::map<int,int> meshPartitioner::getPartToGlobNodeMap(int iPart)
{
  if (iPart > meshParts.size())
  {
    std::cerr << "requested partition number exceeds available partitions" << std::endl;
    exit(1);
  }
  return meshParts[iPart]->getPartToGlobNodeMap();
}

std::map<int,int> meshPartitioner::getPartToGlobElmMap(int iPart)
{
  if (iPart > meshParts.size())
  {
    std::cerr << "requested partition number exceeds available partitions" << std::endl;
    exit(1);
  }
  return meshParts[iPart]->getPartToGlobElmMap();
}
