/* Special purpose class for Proteus HDF5 files */

#include "proteusHdf5.H"
#include <ConversionDriver.H>
#include <TransferBase.H>
#include <TransferDriver.H>
#include <stdio.h>
#include <string.h>
#include <vtkCell.h>
#include <vtkDataSet.h>
#include <vtkFeatureEdges.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkModelMetadata.h>
#include <vtkPointLocator.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkVersion.h>

// MAdLib
#include <MAdLib.h>
#include <NodalDataManager.h>

// GMsh
#include <gmsh/GModel.h>

///////////////////////////////////////////////////
// INITIALIZATION
///////////////////////////////////////////////////

proteusHdf5::proteusHdf5(std::string fname, std::string meshFname,
                         std::string edgeSidesetName, std::string exoMeshFName,
                         bool _lowOrder, bool bndryConst)
    : hdf5Reader(fname, meshFname) {
  // Set boolean determining whether to output higher order elements
  // as lower order elements; useful for visualization
  lowOrder = _lowOrder;

  // Open HDF5 file
  openFile();

  // Get file system metadata
  getControlInformation();
  getVectorNames("ELEMENT_VECTOR_NAMES", elementVectorNames, numElementVectors);
  getVectorNames("VERTEX_VECTOR_NAMES", vertexVectorNames, numVertexVectors);

  // Read in each block
  getBlocks();

  // Merge block data
  mergeBlocks();

  // Set HDF5 data
  setNumberOfDimensions(numDimensions);
  setNumberOfVertices(mySuperBlock.numVertices);
  setNumberOfElements(mySuperBlock.numElements);
  setCoordinates(mySuperBlock.xCrd, mySuperBlock.yCrd, mySuperBlock.zCrd);
  setElementTypesList(mySuperBlock.elementTypesList);
  setElementTypes(mySuperBlock.elementTypes);
  setConnectivities(mySuperBlock.vtkConnectivity);

  // Export to VTK file
  exportToVTKMesh();

  // Export to MeshBase
  exportToMeshBase();

  // Get VTK mesh
  myMeshBase = getMeshBase();

  // Add vertex data to meshBase
  setFields(myMeshBase, vertexVectorNames, mySuperBlock.vertexData, 0);

  // Add element data to meshBase
  setFields(myMeshBase, elementVectorNames, mySuperBlock.elementData, 1);

  // If mesh is 2d, get node ordering convention to follow when splitting cells
  bool rhr;
  // true = right-hand rule (projection of normal vector onto z-axis is +)
  // false = left-hand rule (projection of normal vector onto z-axis is -)
  if (getNumDimensions() == 2) {
    std::cout << "Checking 2d cell node order" << std::endl;
    rhr = get2dCellNodeOrder(myMeshBase);
    if (rhr) {
      std::cout << "2d node ordering right-hand rule: +z" << std::endl;
    } else {
      std::cout << "2d node ordering right-hand rule: -z" << std::endl;
    }
  }

  // Create a meshbase with quads removed (converted to tris)
  std::cout << "Removing quad elements (converting to tris) for refinement"
            << std::endl;
  meshBase *myMeshBaseNoQuads = myMeshBase->convertQuads();

  // Check node ordering again after element conversion (all tris)
  if (getNumDimensions() == 2) {
    std::cout << "Checking 2d cell node order" << std::endl;
    rhr = get2dCellNodeOrder(myMeshBaseNoQuads);
    if (rhr) {
      std::cout << "2d node ordering right-hand rule: +z" << std::endl;
    } else {
      std::cout << "2d node ordering right-hand rule: -z" << std::endl;
    }
  }

  std::cout << "Writing meshBase to file" << std::endl;
  myMeshBase->write();

  // Transfer solution onto new mesh
  std::cout << "Performing solution transfer" << std::endl;
  // myMeshBase->transfer(myMeshBaseNoQuads, "Consistent Interpolation");
  auto transfer = TransferDriver::CreateTransferObject(
      myMeshBase, myMeshBaseNoQuads, "Consistent Interpolation");
  transfer->run();

  // debug
  // myMeshBase->write();
  // myGmshMeshBase->write();

  // Refine mesh
  myMeshBaseNoQuads->refineMesh("uniform", 0.5, "refined.vtu", true,
                                bndryConst);

  // Read in refined mesh
  meshBase *refinedMeshBase = meshBase::Create("refined.vtu");

  // Write refined mesh to file
  myMeshBaseNoQuads->write();
  // myGmshMeshNoQuads->write();

  // Write to file
  // writeVTK();

  // Check node ordering again after refinement
  if (getNumDimensions() == 2) {
    std::cout << "Checking 2d cell node order" << std::endl;
    rhr = get2dCellNodeOrder(refinedMeshBase);
    if (rhr) {
      std::cout << "2d node ordering right-hand rule: +z" << std::endl;
    } else {
      std::cout << "2d node ordering right-hand rule: -z" << std::endl;
    }
  }

  // Sideset implementation

  // Got coordinates for each boundary edge (two pts)
  vtkSmartPointer<vtkPoints> startPts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> endPts = vtkSmartPointer<vtkPoints>::New();
  getBoundaryEdgePts(startPts, endPts);

  // Initialize lists for storing sideset quantities
  vtkSmartPointer<vtkIdList> sidesetElementIdList =
      vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> sidesetSideIdList =
      vtkSmartPointer<vtkIdList>::New();
  // Find elements and sides in meshbase
  // Create two arrays:
  // -one with element ID for sideset cells
  // -one with side ID for each sideset cell
  getBoundarySideSets(refinedMeshBase, startPts, endPts, sidesetElementIdList,
                      sidesetSideIdList);

  // Setting meshbase metadata
  // Set meshbase sideset metadata
  vtkSmartPointer<vtkModelMetadata> sidesetData =
      vtkSmartPointer<vtkModelMetadata>::New();
  vtkSmartPointer<vtkStringArray> sidesetNames =
      vtkSmartPointer<vtkStringArray>::New();
  // currently supports one sideset
  sidesetNames->InsertNextValue(edgeSidesetName);
  sidesetData->SetSideSetNames(sidesetNames);

  // Convert sideset lists into int*
  std::vector<int> _sidesetElementIdList;
  for (int iId = 0; iId < sidesetElementIdList->GetNumberOfIds(); iId++) {
    _sidesetElementIdList.push_back(sidesetElementIdList->GetId(iId));
  }

  std::vector<int> _sidesetSideIdList;
  for (int iId = 0; iId < sidesetSideIdList->GetNumberOfIds(); iId++) {
    _sidesetSideIdList.push_back(sidesetSideIdList->GetId(iId));
  }

  // Set sidesets to meshbase
  sidesetData->SetSideSetElementList(_sidesetElementIdList.data());
  sidesetData->SetSideSetSideList(_sidesetSideIdList.data());
  int sidesetSizes[1] = {
      static_cast<int>(sidesetElementIdList->GetNumberOfIds())};

  sidesetData->SetSideSetSize(sidesetSizes);
  sidesetData->SetNumberOfSideSets(1);
  refinedMeshBase->setMetadata(sidesetData);

  refinedMeshBase->write();

  ConversionDriver *convdrvobj = new ConversionDriver();
  std::vector<meshBase *> mbVec;
  mbVec.push_back(refinedMeshBase);
  std::string exoName = exoMeshFName;

  convdrvobj->genExo(mbVec, exoName);

  // Close HDF5 file
  closeFile();
}

proteusHdf5::~proteusHdf5() {}

void proteusHdf5::getControlInformation() {
  std::cout << "Reading in PROTEUS control data..." << std::endl;

  // Set control buffer information
  H5std_string dsetName("CONTROL");
  int bufSize = 5;
  std::vector<int> buffer(bufSize, 0);

  // Read dataset
  readTopDataSet(dsetName, buffer);

  // Set control data
  numBlocks = buffer[0];
  numDimensions = buffer[1];
  numVertexVectors = buffer[2];
  numElementVectors = buffer[3];
  charStringLength = buffer[4];

  // Output to info to screen
  std::cout << "  Number of Blocks: " << numBlocks << std::endl;
  std::cout << "  Number of Dimensions: " << numDimensions << std::endl;
  std::cout << "  Number of Vertex Vectors: " << numVertexVectors << std::endl;
  std::cout << "  Number of Element Vectors: " << numElementVectors
            << std::endl;
  std::cout << "  Length of Character String: " << charStringLength
            << std::endl;
}

void proteusHdf5::getVectorNames(std::string name,
                                 std::vector<std::string> &stringVector,
                                 int numVectors) {
  std::cout << "Reading in PROTEUS " << name << "..." << std::endl;

  // Set vector buffer information
  H5std_string dsetName(name);
  std::string buffer;

  // Read dataset
  readTopDataSet(dsetName, buffer);

  // Set element vector names
  int bufSize = numVectors * charStringLength;
  std::string tmp = "";
  std::string buf = "";
  int iCharCnt = 0;
  int iChar = 0;
  while (iChar <= bufSize) {
    buf = "";
    buf += buffer[iChar];
    if (!(iCharCnt < charStringLength)) {
      stringVector.push_back(std::string(tmp));
      tmp = "";
      iCharCnt = 0;
    }
    if (buf.compare(" ") != 0) {
      tmp += buffer[iChar];
    }
    iChar++;
    iCharCnt++;
  }

  // Output info to screen
  std::cout << "Found " << name << ":" << std::endl;
  for (auto nameItr = stringVector.begin(); nameItr != stringVector.end();
       nameItr++) {
    std::cout << "  " << *nameItr << std::endl;
  }
}

void proteusHdf5::getBlocks() {
  std::cout << "Reading in PROTEUS block data..." << std::endl;

  for (int iBlock = 0; iBlock < numBlocks; iBlock++) {
    // Set PROTEUS convention block name
    std::string blockName = "BLOCK";
    std::string numStr = std::to_string(iBlock + 1); // Proteus blocks 1-indexed
    while (numStr.length() < 12) {
      numStr.insert(0, "0");
    }
    blockName += numStr;
    std::cout << "Reading in " + blockName + "..." << std::endl;

    // Create HDF5 Group object for block
    Group group;
    readGroup(blockName, group);

    // Create block object
    proteusBlock myBlock = proteusBlock();
    myBlock.blockName = blockName;

    // Read block info
    getBlockInfo(group, myBlock);

    // Read coordinates
    getBlockXYZ(group, myBlock);

    // Read global IDs
    getBlockGlobalID(group, myBlock);

    // Read element data
    getBlockElementData(group, myBlock);

    // Read vertex data
    getBlockVertexData(group, myBlock);

    // Push back proteus block into data structure
    proteusBlocks.push_back(myBlock);
  }

  std::cout << "Finished reading PROTEUS block data..." << std::endl;
}

void proteusHdf5::mergeBlocks() {
  std::cout << "Merging PROTEUS block data..." << std::endl;

  // Determine total number of vertices based on maximum global ID
  int numGlobalVertices = 0;

  // Determine total number of elements
  int numGlobalElements = 0;
  int maxNumVerticesPerElement = 0;

  // Determine number of element types
  std::vector<int> elementTypesList;
  std::vector<int> originalElementTypesList;

  for (auto blockItr = proteusBlocks.begin(); blockItr != proteusBlocks.end();
       blockItr++) {
    // Get number of global vertex Ids
    int localMaxVertexId = *std::max_element((*blockItr).loc2Glob.begin(),
                                             (*blockItr).loc2Glob.end());
    if (localMaxVertexId > numGlobalVertices) {
      numGlobalVertices = localMaxVertexId;
    }

    // Accumulate number of global element Ids
    numGlobalElements += (*blockItr).numElements;

    // Get maximum number of vertices per element
    if ((*blockItr).numVerticesPerElement > maxNumVerticesPerElement) {
      maxNumVerticesPerElement = (*blockItr).numVerticesPerElement;
    }

    if (std::find(originalElementTypesList.begin(),
                  originalElementTypesList.end(),
                  (*blockItr).originalVtkElementType) ==
        originalElementTypesList.end()) {
      originalElementTypesList.push_back((*blockItr).originalVtkElementType);
      elementTypesList.push_back((*blockItr).vtkElementType);
    }
  }
  numGlobalVertices++; // global vertices are 0-indexed, so total number is
                       // maxID+1

  // Set superblock info
  mySuperBlock.numVertices = numGlobalVertices;
  mySuperBlock.numElements = numGlobalElements;
  mySuperBlock.maxNumVerticesPerElement = maxNumVerticesPerElement;

  std::cout << "    Found " << numGlobalVertices << " global vertices"
            << std::endl;
  std::cout << "    Found " << numGlobalElements << " global elements"
            << std::endl;

  // Vector to store global vertices
  std::vector<std::vector<double>> coordinates(
      numGlobalVertices, std::vector<double>(numDimensions, 0.0));

  // Vector to store global elements
  std::vector<std::vector<int>> elements(
      numGlobalElements, std::vector<int>(maxNumVerticesPerElement, 0));

  // Map of element to VTK element type
  std::vector<int> elementTypes(numGlobalElements);

  // Vector to store vertex data
  std::vector<std::vector<double>> vertexData(
      numVertexVectors, std::vector<double>(numGlobalVertices, 0.0));

  // Vector to store element data
  std::vector<std::vector<double>> elementData(
      numElementVectors, std::vector<double>(numGlobalElements, 0.0));

  // Instantiate variables for loops below
  int locVertId, globElemId, locElemId;
  int lId, gId, vertNo, dimNo, origId;
  globElemId = 0;

  // Loop over blocks
  std::cout << "    Accumulating mesh information from each block..."
            << std::endl;
  for (auto blockItr = proteusBlocks.begin(); blockItr != proteusBlocks.end();
       blockItr++) {
    std::cout << "        " << (*blockItr).blockName << ":" << std::endl;
    locVertId = 0;
    locElemId = 0;
    // Loop over element coordinates
    for (auto elemItr = (*blockItr).vertices.begin();
         elemItr != (*blockItr).vertices.end(); elemItr++) {
      // Loop over dimension coordinates
      origId = locVertId;
      for (auto dimItr = (*elemItr).begin(); dimItr != (*elemItr).end();
           dimItr++) {
        locVertId = origId;
        dimNo = distance((*elemItr).begin(), dimItr);
        // Loop over vertex coordinates
        for (auto vertItr = (*dimItr).begin(); vertItr != (*dimItr).end();
             vertItr++) {
          vertNo = distance((*dimItr).begin(), vertItr);
          elements[globElemId][vertNo] = (*blockItr).loc2Glob[locVertId];
          coordinates[(*blockItr).loc2Glob[locVertId]][dimNo] = *vertItr;
          locVertId++;
        }
        elementTypes[globElemId] = (*blockItr).vtkElementType;
      }

      // Accumulate element field data
      for (int iElemData = 0; iElemData < numElementVectors; iElemData++) {
        elementData[iElemData][globElemId] =
            (*blockItr).elementData[iElemData][locElemId];
      }

      locElemId++;
      globElemId++;
    }

    // Loop over vertices
    for (auto vertItr = (*blockItr).loc2Glob.begin();
         vertItr != (*blockItr).loc2Glob.end(); vertItr++) {
      lId = distance((*blockItr).loc2Glob.begin(), vertItr);
      gId = *vertItr;
      for (int iVertData = 0; iVertData < numVertexVectors; iVertData++) {
        vertexData[iVertData][gId] = (*blockItr).vertexData[iVertData][lId];
      }
    }

    std::cout << "            " << (*blockItr).vertices.size() << " vertices"
              << std::endl;
    std::cout << "            " << (*blockItr).numElements << " elements"
              << std::endl;
    std::cout << "            " << numVertexVectors << " vertex fields"
              << std::endl;
    std::cout << "            " << numElementVectors << " element fields"
              << std::endl;
  }

  // Re-structure data into separate x,y,z coordinates (todo: could be done in
  // loops above)
  std::cout << "    Separating coordinate arrays into individual dimension "
               "arrays..."
            << std::endl;
  for (auto vertItr = coordinates.begin(); vertItr != coordinates.end();
       vertItr++) {
    mySuperBlock.xCrd.push_back((*vertItr)[0]);
    if (numDimensions > 1) {
      mySuperBlock.yCrd.push_back((*vertItr)[1]);
    }
    if (numDimensions > 2) {
      mySuperBlock.zCrd.push_back((*vertItr)[2]);
    }
  }

  // Convert element connectivity information into a format that is easier to
  // parse for hdf5Reader class. (todo: chould be done in loops above) Main
  // structure shown below: elementConnectivity[element type no][global
  // element id][vertex no] = [global vertex id] elementTypesList[element type
  // no] = vtk element type integer
  std::cout << "    Sorting element connectivities by element type..."
            << std::endl;
  std::vector<std::map<int, std::vector<int>>> elementConnectivity;

  // Loop over original element type list (contains all element types read in
  // across all blocks)
  auto origTypeItr = originalElementTypesList.begin();
  while (origTypeItr != originalElementTypesList.end()) {
    // temporary map for storing connectivity for each element type
    std::map<int, std::vector<int>> tmpElementConnectivity;

    int vertNum, elemId;
    // Loop over elements and their types
    auto elemItr = elements.begin();
    auto elemTypeItr = elementTypes.begin();
    while (elemItr != elements.end() || elemTypeItr != elementTypes.end()) {
      // Get global element id
      elemId = distance(elements.begin(), elemItr);
      // Loop over all vertices of each element
      for (auto vertItr = (*elemItr).begin(); vertItr != (*elemItr).end();
           vertItr++) {
        // Get vertex number
        vertNum = distance((*elemItr).begin(), vertItr);
        // Check element type against each supported type
        if (*elemTypeItr == VTK_QUAD) {
// Remove extra vertices from high-order Lagrange quads if output
// to lower-order quad is desired
#if VTK_MAJOR_VERSION > 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 1)
          if (*origTypeItr == VTK_LAGRANGE_QUADRILATERAL) {
            if (vertNum == 0 || vertNum == 3 || vertNum == 12 ||
                vertNum == 15) {
              tmpElementConnectivity[elemId].push_back(*vertItr);
            }
          } else {
            tmpElementConnectivity[elemId].push_back(*vertItr);
          }
#else
          if (*origTypeItr == -1) // lower VTK versions don't have
                                  // LAGRANGE_QUADRILATERAL types
          {
            if (vertNum == 0 || vertNum == 3 || vertNum == 12 ||
                vertNum == 15) {
              tmpElementConnectivity[elemId].push_back(*vertItr);
            }
          } else {
            tmpElementConnectivity[elemId].push_back(*vertItr);
          }
#endif
        }
#if VTK_MAJOR_VERSION > 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 1)
        else if (*elemTypeItr == VTK_LAGRANGE_QUADRILATERAL) {
          tmpElementConnectivity[elemId].push_back(*vertItr);
        }
#endif
        else if (*elemTypeItr == VTK_TRIANGLE) {
          tmpElementConnectivity[elemId].push_back(*vertItr);
        }
      }

      // Increment element and element type iterators
      if (elemItr != elements.end()) {
        ++elemItr;
      }
      if (elemTypeItr != elementTypes.end()) {
        ++elemTypeItr;
      }
    }

    // Push back each inidividual element type connectivity
    // into vector containing all types
    elementConnectivity.push_back(tmpElementConnectivity);

    if (origTypeItr != originalElementTypesList.end()) {
      ++origTypeItr;
    }
  }

  // Re-order vertices per VTK convention
  std::cout << "    Re-ordering element vertices per VTK convention..."
            << std::endl;
  int vertCount;
  int swapVert1 = -1, swapVert2 = -1;
  int connId;
  // If VTK_QUAD, re-order vertices
  auto testIt =
      std::find(elementTypesList.begin(), elementTypesList.end(), VTK_QUAD);
  if (testIt != elementTypesList.end()) {
    std::cout << "        Re-ordering VTK_QUAD element vertices..."
              << std::endl;
    connId = distance(elementTypesList.begin(), testIt);
    for (auto elemIt = elementConnectivity[connId].begin();
         elemIt != elementConnectivity[connId].end(); elemIt++) {
      vertCount = 0;
      for (int i = 0; i < (elemIt->second).size(); i++) {
        if (vertCount == 2) {
          swapVert1 = (elemIt->second)[i];
        }
        if (vertCount == 3) {
          swapVert2 = (elemIt->second)[i];
          (elemIt->second)[i] = swapVert1;
          (elemIt->second)[i - 1] = swapVert2;
          vertCount = -1;
        }
        vertCount++;
      }
    }
  }
  // If VTK_TRI, re-order vertices
  /*
  testIt = std::find(elementTypesList.begin(), elementTypesList.end(),
  VTK_TRIANGLE); if (testIt != elementTypesList.end())
  {
    std::cout << "        Re-ordering VTK_TRIANGLE element vertices..." <<
  std::endl; connId = distance(elementTypesList.begin(), testIt); for (auto
  elemIt = elementConnectivity[connId].begin(); elemIt !=
  elementConnectivity[connId].end(); elemIt++)
    {
      vertCount = 0;
      for (int i = 0; i < (elemIt->second).size(); i++)
      {
        if (vertCount == 1)
        {
          swapVert1 = (elemIt->second)[i];
        }
        if (vertCount == 2)
        {
          swapVert2 = (elemIt->second)[i];
          (elemIt->second)[i] = swapVert1;
          (elemIt->second)[i-1] = swapVert2;
          vertCount = -1;
        }
        vertCount++;
      }
    }
  }
  */

  // Set data in superBlock struct
  std::cout << "    Setting mesh data in PROTEUS superblock..." << std::endl;
  mySuperBlock.coordinates = coordinates;
  mySuperBlock.elements = elements;
  mySuperBlock.vertexData = vertexData;
  mySuperBlock.elementData = elementData;
  mySuperBlock.vtkConnectivity = elementConnectivity;
  mySuperBlock.elementTypesList = elementTypesList;
  mySuperBlock.originalElementTypesList = originalElementTypesList;
  mySuperBlock.elementTypes = elementTypes;

  std::cout << "Finished merging PROTEUS block data..." << std::endl;
}

void proteusHdf5::getBlockInfo(Group &group, proteusBlock &myBlock) {
  int bufSize = 3;
  std::vector<int> buffer(bufSize, 0);
  readGroupDataSet("INFO", group, buffer);

  std::vector<int> tmpVec;
  for (int iBuf = 0; iBuf < bufSize; iBuf++) {
    tmpVec.push_back(buffer[iBuf]);
  }
  myBlock.numElements = tmpVec[0];
  myBlock.numVerticesPerElement = tmpVec[1];

  // parse element type
  // Note: the element ID numbers below don't match up with the PROTEUS
  // documentation, but have been implemented based on PROTEUS HDF5 test files
  // provided
  if (tmpVec[2] == 162) {
    if (lowOrder) {
#if VTK_MAJOR_VERSION > 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 1)
      myBlock.originalVtkElementType = VTK_LAGRANGE_QUADRILATERAL;
#else
      myBlock.originalVtkElementType = -1; // older VTK versions don't support
                                           // LAGRANGE_QUADRILATERAL types
#endif
      myBlock.vtkElementType = VTK_QUAD;
    } else {
#if VTK_MAJOR_VERSION > 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 1)
      myBlock.originalVtkElementType = VTK_LAGRANGE_QUADRILATERAL;
      myBlock.vtkElementType = VTK_LAGRANGE_QUADRILATERAL;
#else
      myBlock.originalVtkElementType = VTK_QUAD;
      myBlock.vtkElementType = VTK_QUAD;
      std::cout << "Warning: VTK version " << vtkVersion::GetVTKSourceVersion()
                << " does not support Lagrange Quadrilateral cells."
                << std::endl;
#endif
    }
  } else if (tmpVec[2] == 160) {
    myBlock.originalVtkElementType = VTK_QUAD;
    myBlock.vtkElementType = VTK_QUAD;
  }
  // else if (tmpVec[2] == 110)
  // It looks like some of the triangle elements in Proteus test files have an
  // element type of 101 instead of 110, as documented in the user manual. I
  // have just added this as an additional element type to catch this.
  else if (tmpVec[2] == 110) {
    myBlock.originalVtkElementType = VTK_TRIANGLE;
    myBlock.vtkElementType = VTK_TRIANGLE;
  } else {
    std::cerr << "Element type not supported: " << tmpVec[2] << std::endl;
    exit(-1);
  }

  // Output info to screen
  std::cout << "   Info:" << std::endl;
  std::cout << "      Number of elements: " << tmpVec[0] << std::endl;
  std::cout << "      Number of vertices per element: "
            << myBlock.numVerticesPerElement << std::endl;
  std::cout << "      VTK element type: " << myBlock.vtkElementType
            << std::endl;
}

void proteusHdf5::getBlockXYZ(Group &group, proteusBlock &myBlock) {
  // Define xyz buffer
  int bufSize =
      numDimensions * myBlock.numElements * myBlock.numVerticesPerElement;
  std::vector<float> buffer(bufSize, 0);

  // Read dataset from group
  readGroupDataSet("XYZ", group, buffer);

  // Output info to screen
  std::cout << "   Coordinate info:" << std::endl;
  std::cout << "      Read in " << bufSize << " vertex coordinates"
            << std::endl;
  for (int iDim = 0; iDim < numDimensions; iDim++) {
    std::cout << "      " << myBlock.numElements * myBlock.numVerticesPerElement
              << " in dimension " << iDim << std::endl;
  }

  // Convert vertex coordinates into an un-flattened array using row-major
  // ordering
  std::vector<std::vector<std::vector<double>>> tmpVec(
      myBlock.numElements,
      std::vector<std::vector<double>>(
          numDimensions,
          std::vector<double>(myBlock.numVerticesPerElement, 0.0)));
  int tmpInd = 0;
  for (int iElem = 0; iElem < myBlock.numElements; iElem++) {
    for (int iDim = 0; iDim < numDimensions; iDim++) {
      for (int iVert = 0; iVert < myBlock.numVerticesPerElement; iVert++) {
        tmpVec[iElem][iDim][iVert] = double(buffer[tmpInd]);
        tmpInd++;
      }
    }
  }

  // Assign vertex coordinates to block
  myBlock.vertices = tmpVec;
}

void proteusHdf5::getBlockGlobalID(Group &group, proteusBlock &myBlock) {
  // Define GlobalID buffer
  int bufSize = myBlock.numElements * myBlock.numVerticesPerElement;
  std::vector<int> buffer(bufSize, 0);

  // Read dataset from group
  readGroupDataSet("GLOBALID", group, buffer);

  // Output info to screen
  std::cout << "   Global ID info:" << std::endl;
  std::cout << "      Read in "
            << myBlock.numElements * myBlock.numVerticesPerElement
            << " global IDs" << std::endl;

  // Create a map from local to global ID
  std::vector<int> tmpVec;
  for (int iBuf = 0; iBuf < bufSize; iBuf++) {
    tmpVec.push_back(buffer[iBuf]);
  }

  // Assign local to global ID map to block
  myBlock.loc2Glob = tmpVec;
}

void proteusHdf5::getBlockElementData(Group &group, proteusBlock &myBlock) {
  // Define GlobalID buffer
  int bufDim1 = numElementVectors;
  int bufDim2 = myBlock.numElements;
  int bufDimFlat = bufDim1 * bufDim2;

  // Buffer for HDF read (float)
  std::vector<float> flatBuffer(bufDimFlat, 0.0);

  // Buffer for proteusHDF5 object (2d double)
  std::vector<std::vector<double>> buffer(bufDim1,
                                          std::vector<double>(bufDim2, 0.0));

  // Read dataset from group
  readGroupDataSet("ELEMENTDATA", group, flatBuffer);

  // Convert to double
  std::vector<double> doubleFlatBuffer(flatBuffer.begin(), flatBuffer.end());

  // Unflatten buffer to 2d
  unflattenBuffer2d(doubleFlatBuffer, buffer);
  // Todo: Make this a generic templated function like flattenBuffer() above
  // auto resizedBuffer = unflattenBuffer(flatBuffer);

  // Output info to screen
  std::cout << "   Element Data info:" << std::endl;
  std::cout << "      Read in " << numElementVectors
            << " element data:" << std::endl;
  for (auto edItr = elementVectorNames.begin();
       edItr != elementVectorNames.end(); edItr++) {
    std::cout << "      " << *edItr << std::endl;
  }

  // Assign element data to block
  myBlock.elementData = buffer;
}

void proteusHdf5::getBlockVertexData(Group &group, proteusBlock &myBlock) {
  // Define GlobalID buffer
  int bufDim1 = numVertexVectors;
  int bufDim2 = myBlock.numElements * myBlock.numVerticesPerElement;
  int bufDimFlat = bufDim1 * bufDim2;

  // Buffer for HDF read (float)
  std::vector<float> flatBuffer(bufDimFlat, 0.0);

  // Buffer for proteusHDF5 object (2d double)
  std::vector<std::vector<double>> buffer(bufDim1,
                                          std::vector<double>(bufDim2, 0.0));

  // Read dataset from group
  readGroupDataSet("VERTEXDATA", group, flatBuffer);

  // Convert to double
  std::vector<double> doubleFlatBuffer(flatBuffer.begin(), flatBuffer.end());

  // Unflatten buffer to 2d
  unflattenBuffer2d(doubleFlatBuffer, buffer);

  // Output info to screen
  std::cout << "   Vertex Data info:" << std::endl;
  std::cout << "      Read in " << numVertexVectors
            << " vertex data:" << std::endl;
  for (auto vdItr = vertexVectorNames.begin(); vdItr != vertexVectorNames.end();
       vdItr++) {
    std::cout << "      " << *vdItr << std::endl;
  }

  // Assign vertex data to block
  myBlock.vertexData = buffer;
}

// unused function for writing boundary condition information to json file
/*
void proteusHdf5::writeBcFile(meshBase* myMeshBase, std::string
edgeSidesetName)
{
  if (numDimensions != 2)
  {
    std::cerr << "Error: Sidesets only supports for 2D meshes" << std::endl;
    exit(-1);
  }

  // Construct json data
  json outJson;
  outJson["Operation"] = "Boundary Condition Assignment";
  outJson["Input Mesh"] = "refined.msh";

  json condJson;
  condJson["Name"] = "side_set_0000001";
  condJson["Boundary Type"] = "Edges";
  condJson["Condition Type"] = "Fixed";

  json startArray = json::array();
  json endArray = json::array();
  json startSubArray = json::array();
  json endSubArray = json::array();

  MAd::pGModel tmpMdl = NULL;
  MAd::GM_create(&tmpMdl,"");
  MAd::pMesh refinedMadMesh = MAd::M_new(tmpMdl);
  MAd::M_load(refinedMadMesh,"refined.msh");

  MAd::GM_create(&tmpMdl,"");
  MAd::pMesh skinnedMadMesh = MAd::M_new(tmpMdl);

  std::vector<int> regIds;
  refinedMadMesh->skin_me(skinnedMadMesh, regIds);

  MAd::EIter eit = MAd::M_edgeIter(skinnedMadMesh);
  MAd::pEdge pe;
  while ((pe = MAd::EIter_next(eit)))
  {
    MAd::pVertex v = MAd::E_vertex(pe, 0);
    startSubArray = json::array();
    startSubArray.add(v->X);
    startSubArray.add(v->Y);
    startSubArray.add(v->Z);
    startArray.add(startSubArray);

    v = MAd::E_vertex(pe, 1);
    endSubArray = json::array();
    endSubArray.add(v->X);
    endSubArray.add(v->Y);
    endSubArray.add(v->Z);
    endArray.add(endSubArray);
  }

  condJson["Params"]["Start"] = startArray;
  condJson["Params"]["End"] = endArray;

  json condArray = json::array();
  condArray.add(condJson);

  outJson["Condition"] = condArray;

  ofstream of("./bcs.json");
  of << pretty_print(outJson);
  of.close();
}
*/

void proteusHdf5::getBoundaryEdgePts(vtkSmartPointer<vtkPoints> startPts,
                                     vtkSmartPointer<vtkPoints> endPts) {
  if (numDimensions != 2) {
    std::cerr << "Error: Sidesets only supports for 2D meshes" << std::endl;
    exit(-1);
  }

  MAd::pGModel tmpMdl = NULL;
  MAd::GM_create(&tmpMdl, "");
  MAd::pMesh refinedMadMesh = MAd::M_new(tmpMdl);
  MAd::M_load(refinedMadMesh, "refined.msh");

  MAd::GM_create(&tmpMdl, "");
  MAd::pMesh skinnedMadMesh = MAd::M_new(tmpMdl);

  std::vector<int> regIds;
  refinedMadMesh->skin_me(skinnedMadMesh, regIds);

  MAd::EIter eit = MAd::M_edgeIter(skinnedMadMesh);
  MAd::pEdge pe;
  int edgeId = 0;
  double coord[3];
  while ((pe = MAd::EIter_next(eit))) {
    MAd::pVertex v = MAd::E_vertex(pe, 0);
    coord[0] = v->X;
    coord[1] = v->Y;
    coord[2] = v->Z;
    startPts->InsertPoint(edgeId, coord);

    v = MAd::E_vertex(pe, 1);
    coord[0] = v->X;
    coord[1] = v->Y;
    coord[2] = v->Z;
    endPts->InsertPoint(edgeId, coord);

    edgeId++;
  }
}

void proteusHdf5::getBoundarySideSets(
    meshBase *myMeshBase, vtkSmartPointer<vtkPoints> startPts,
    vtkSmartPointer<vtkPoints> endPts,
    vtkSmartPointer<vtkIdList> sidesetElementIdList,
    vtkSmartPointer<vtkIdList> sidesetSideIdList) {
  // Build up locator
  vtkSmartPointer<vtkPointLocator> pointLocator =
      vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet(myMeshBase->getDataSet());
  pointLocator->BuildLocator();

  // edge Pt iterators
  vtkIdType id1, id2;
  vtkSmartPointer<vtkIdList> cellList1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellList2 = vtkSmartPointer<vtkIdList>::New();
  vtkCell *foundCell;
  vtkSmartPointer<vtkIdList> edgeIdListCompare =
      vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> edgeIdListOriginal =
      vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> edgeIdListIntersection =
      vtkSmartPointer<vtkIdList>::New();
  for (int iBndEdge = 0; iBndEdge < startPts->GetNumberOfPoints(); iBndEdge++) {
    // std::cout << "boundary edge: " << iBndEdge << std::endl;
    edgeIdListOriginal->Reset();

    id1 = pointLocator->FindClosestPoint(startPts->GetPoint(iBndEdge));
    id2 = pointLocator->FindClosestPoint(endPts->GetPoint(iBndEdge));

    edgeIdListOriginal->InsertNextId(id1);
    edgeIdListOriginal->InsertNextId(id2);

    myMeshBase->getDataSet()->GetPointCells(id1, cellList1);
    myMeshBase->getDataSet()->GetPointCells(id2, cellList2);

    cellList1->IntersectWith(cellList2);
    if (cellList1->GetNumberOfIds() != 1) {
      std::cerr << "Error: found multiple cells for boundary edge in a 2d mesh"
                << std::endl;
      exit(-1);
    } else {
      sidesetElementIdList->InsertNextId(cellList1->GetId(0));
      // std::cout << "    elem ID = " << cellList1->GetId(0) << std::endl;
      foundCell = myMeshBase->getDataSet()->GetCell(cellList1->GetId(0));
      vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
      myMeshBase->getDataSet()->GetCellPoints(cellList1->GetId(0), cellPts);
      for (int ipt = 0; ipt < cellPts->GetNumberOfIds(); ipt++) {
        double position[3];
        myMeshBase->getDataSet()->GetPoint(cellPts->GetId(ipt), position);
        // std::cout << "    ipt " << ipt << " position = " << position[0] << ",
        // "
        //          << position[1] << ", " << position[2] << std::endl;
      }
    }

    for (int iCellEdge = 0; iCellEdge < foundCell->GetNumberOfEdges();
         iCellEdge++) {
      edgeIdListCompare = foundCell->GetEdge(iCellEdge)->GetPointIds();

      // Alternate 'intersection' implementation due to pointer issues
      edgeIdListIntersection->Reset();
      for (int iId = 0; iId < edgeIdListOriginal->GetNumberOfIds(); iId++) {
        for (int jId = 0; jId < edgeIdListCompare->GetNumberOfIds(); jId++) {
          if (edgeIdListOriginal->GetId(iId) == edgeIdListCompare->GetId(jId)) {
            edgeIdListIntersection->InsertNextId(
                edgeIdListOriginal->GetId(iId));
          }
        }
      }

      if (edgeIdListIntersection->GetNumberOfIds() == 2) {
        sidesetSideIdList->InsertNextId(iCellEdge);
        // std::cout << "    side ID = " << iCellEdge << std::endl;
        break;
      }
      if (iCellEdge == foundCell->GetNumberOfEdges() - 1) {
        std::cout << "Couldn't find a matching side!" << std::endl;
      }
    }
  }
}

int proteusHdf5::getNumBlocks() const { return (numBlocks); }

int proteusHdf5::getNumDimensions() const { return (numDimensions); }

int proteusHdf5::getNumVertexVectors() const { return (numVertexVectors); }

int proteusHdf5::getNumElementVectors() const { return (numElementVectors); }

int proteusHdf5::getCharStringLength() const { return (charStringLength); }

bool proteusHdf5::get2dCellNodeOrder(meshBase *myMeshBase) {
  std::map<int, bool> rhr;
  for (int iCell = 0; iCell < myMeshBase->getDataSet()->GetNumberOfCells();
       iCell++) {
    // Get first cell
    vtkSmartPointer<vtkCell> testCell =
        myMeshBase->getDataSet()->GetCell(iCell);
    vtkSmartPointer<vtkPoints> cellPts = testCell->GetPoints();
    double pt0[3];
    double pt1[3];
    double pt2[3];
    double vc0[3];
    double vc1[3];
    double vc2[3];
    double zNorm[3] = {0, 0, 1};
    double res;
    if (testCell->GetCellType() == VTK_TRIANGLE) {
      cellPts->GetPoint(0, pt0);
      cellPts->GetPoint(1, pt1);
      cellPts->GetPoint(2, pt2);
      vtkMath::Subtract(pt1, pt0, vc0);
      vtkMath::Subtract(pt2, pt0, vc1);
    } else if (testCell->GetCellType() == VTK_QUAD) {
      cellPts->GetPoint(0, pt0);
      cellPts->GetPoint(1, pt1);
      cellPts->GetPoint(2, pt2);
      vtkMath::Subtract(pt1, pt0, vc0);
      vtkMath::Subtract(pt2, pt0, vc1);
    } else if (testCell->GetCellType() != VTK_QUAD &&
               testCell->GetCellType() != VTK_TRIANGLE) {
      std::cerr << "Only triangular and quadrilateral elements supported."
                << std::endl;
      exit(-1);
    } else {
      continue;
    }
    vtkMath::Cross(vc0, vc1, vc2);
    res = vtkMath::Dot(vc2, zNorm);

    bool ans;
    if (res > 0) {
      ans = true;
    } else {
      ans = false;
    }
    if ((testCell->GetCellType() == VTK_QUAD &&
         rhr.find(VTK_QUAD) == rhr.end()) ||
        (testCell->GetCellType() == VTK_TRIANGLE &&
         rhr.find(VTK_TRIANGLE) == rhr.end())) {
      rhr[testCell->GetCellType()] = ans;
    } else {
      if (rhr[testCell->GetCellType()] != ans) {
        std::cout << "Mismatched ordering for type " << testCell->GetCellType()
                  << ", cell number " << iCell << std::endl;
      }
    }
  }
  bool conv = false;
  for (auto itr = rhr.begin(); itr != rhr.end(); itr++) {
    if (itr == rhr.begin()) {
      conv = itr->second;
    } else {
      if (conv != itr->second) {
        std::cerr << "Error: Mesh contains multiple element types with "
                     "different node ordering."
                  << std::endl;
        exit(-1);
      }
    }
  }
  return conv;
}
