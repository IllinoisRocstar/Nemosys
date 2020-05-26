#include <hdf5Reader.H>
#include <meshBase.H>
#include <string.h>
#include <algorithm>
#include <iostream>

/********************************************
    hdfReader class implementation
*********************************************/

std::string hdf5Reader::getFileName() { return h5FileName; }

int hdf5Reader::openFile() {
  std::cout << "Opening HDF5 file " << h5FileName << std::endl;

  try {
    // Open HDF5 file
    H5std_string h5FileName_(h5FileName);
    file.openFile(h5FileName_, H5F_ACC_RDWR);
  } catch (FileIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 11 && H5_VERS_RELEASE >= 0
      error.printErrorStack();
#else
    error.printError();
#endif
    return -1;
  }
  return 0;
}

int hdf5Reader::closeFile() {
  std::cout << "Closing HDF5 file " << h5FileName << std::endl;

  try {
    // Close HDF5 file
    file.close();
  } catch (FileIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 11 && H5_VERS_RELEASE >= 0
      error.printErrorStack();
#else
    error.printError();
#endif
    return -1;
  }
  return 0;
}

// Read string data
int hdf5Reader::readTopDataSet(std::string dsetName, std::string &buffer) {
  try {
    // Open dataset
    DataSet dset = file.openDataSet(dsetName);
    // Read dataset
    readDataSet(dset, buffer);
  }
  // catch failure caused by the H5File operations
  catch (FileIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 11 && H5_VERS_RELEASE >= 0
      error.printErrorStack();
#else
    error.printError();
#endif
    return -1;
  }
  return 0;
}

// Read string data
int hdf5Reader::readDataSet(DataSet dset, std::string &buffer) {
  try {
    // Get type of data
    DataType dtype = dset.getDataType();
    // Read data
    dset.read(buffer, dtype);
    // Close dataset
    dset.close();
  }
  // catch failure caused by the H5File operations
  catch (DataSetIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 11 && H5_VERS_RELEASE >= 0
      error.printErrorStack();
#else
    error.printError();
#endif
    return -1;
  }
  return 0;
}

// Read existing HDF5 Group
int hdf5Reader::readGroup(std::string groupName, Group &group) {
  try {
    // Open group
    group = Group(file.openGroup(groupName));
  } catch (FileIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
    error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 11 && H5_VERS_RELEASE >= 0
      error.printErrorStack();
#else
    error.printError();
#endif
    return -1;
  }
  return 0;
}

void hdf5Reader::setNumberOfVertices(int _numVertices) {
  numVertices = _numVertices;
}

void hdf5Reader::setNumberOfElements(int _numElements) {
  numElements = _numElements;
}

void hdf5Reader::setNumberOfDimensions(int _numDimensions) {
  numDimensions = _numDimensions;
}

void hdf5Reader::setCoordinates(std::vector<double> _xCrd) { xCrd = _xCrd; }

void hdf5Reader::setCoordinates(std::vector<double> _xCrd,
                                std::vector<double> _yCrd) {
  xCrd = _xCrd;
  yCrd = _yCrd;
}

void hdf5Reader::setCoordinates(std::vector<double> _xCrd,
                                std::vector<double> _yCrd,
                                std::vector<double> _zCrd) {
  xCrd = _xCrd;
  yCrd = _yCrd;
  zCrd = _zCrd;
}

void hdf5Reader::setElementTypes(std::vector<int> _vtkElementTypes) {
  vtkElementTypes = _vtkElementTypes;
}

void hdf5Reader::setElementTypesList(std::vector<int> _vtkElementTypesList) {
  vtkElementTypesList = _vtkElementTypesList;
}

void hdf5Reader::setConnectivities(
    std::vector<std::map<int, std::vector<int>>> _elementConnectivity) {
  elementConnectivity = _elementConnectivity;
}

std::vector<int> hdf5Reader::getElementConnectivity(int elementId) {
  int nVerts;
  int connId = -1;
#if VTK_MAJOR_VERSION > 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 1)
  if (vtkElementTypes[elementId] == VTK_LAGRANGE_QUADRILATERAL) {
    nVerts = 16;
    auto itr = std::find(vtkElementTypesList.begin(), vtkElementTypesList.end(),
                         VTK_LAGRANGE_QUADRILATERAL);
    if (itr != vtkElementTypesList.end()) {
      connId = distance(vtkElementTypesList.begin(), itr);
    }
  } else if (vtkElementTypes[elementId] == VTK_QUAD)
#else
  if (vtkElementTypes[elementId] == VTK_QUAD)
#endif
  {
    nVerts = 4;
    auto itr = std::find(vtkElementTypesList.begin(), vtkElementTypesList.end(),
                         VTK_QUAD);
    if (itr != vtkElementTypesList.end()) {
      connId = distance(vtkElementTypesList.begin(), itr);
    }
  } else if (vtkElementTypes[elementId] == VTK_TRIANGLE) {
    nVerts = 3;
    auto itr = std::find(vtkElementTypesList.begin(), vtkElementTypesList.end(),
                         VTK_TRIANGLE);
    if (itr != vtkElementTypesList.end()) {
      connId = distance(vtkElementTypesList.begin(), itr);
    }
  } else {
    std::cerr << "Element type not supported: " << vtkElementTypes[elementId]
              << std::endl;
    exit(-1);
  }
  if (connId != -1) {
    std::vector<int> subConnectivity(
        elementConnectivity[connId][elementId].begin(),
        elementConnectivity[connId][elementId].begin() + nVerts);
    return subConnectivity;
  } else {
    std::cerr << "Element type " << vtkElementTypes[elementId]
              << " not found in connectivity table" << std::endl;
    exit(-1);
  }
}

void hdf5Reader::exportToVTKMesh() {
  std::cout << "Exporting hdf5Reader object to VTK format..." << std::endl;

  if (!vtkMesh) {
    // declare vtk dataset
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    // points to be pushed into dataSet
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // set double precision points
    points->SetDataTypeToDouble();

    // allocate size for vtk point container
    points->SetNumberOfPoints(numVertices);

    if (numDimensions == 2) {
      for (int i = 0; i < numVertices; ++i) {
        points->SetPoint(i, xCrd[i], yCrd[i], 0.0);
      }
    } else {
      for (int i = 0; i < numVertices; ++i) {
        points->SetPoint(i, xCrd[i], yCrd[i], zCrd[i]);
      }
    }
    // add points to vtk mesh data structure
    std::cout << "    setting VTK point data..." << std::endl;
    dataSet_tmp->SetPoints(points);
    // allocate space for elements
    std::cout << "    allocating memory for elements..." << std::endl;
    dataSet_tmp->Allocate(numElements);
    // add the elements
    std::cout << "    adding VTK elements..." << std::endl;
    for (int i = 0; i < numElements; ++i) {
      vtkSmartPointer<vtkIdList> vtkElmIds = vtkSmartPointer<vtkIdList>::New();
      std::vector<int> elmIds(getElementConnectivity(i));
      vtkElmIds->SetNumberOfIds(elmIds.size());
      for (int j = 0; j < elmIds.size(); ++j) {
        vtkElmIds->SetId(j, elmIds[j] - 0);
      }
      switch (vtkElementTypes[i]) {
#if VTK_MAJOR_VERSION > 8 || (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION >= 1)
        case VTK_LAGRANGE_QUADRILATERAL:
          dataSet_tmp->InsertNextCell(VTK_LAGRANGE_QUADRILATERAL, vtkElmIds);
          break;
#endif
        case VTK_QUAD:
          dataSet_tmp->InsertNextCell(VTK_QUAD, vtkElmIds);
          break;
        case VTK_TRIANGLE:
          dataSet_tmp->InsertNextCell(VTK_TRIANGLE, vtkElmIds);
          break;
        default:
          std::cerr << "Unknown element type " << vtkElementTypes[i]
                    << std::endl;
          break;
      }
    }
    vtkMesh = dataSet_tmp;
  }
}

vtkSmartPointer<vtkDataSet> hdf5Reader::getVTKMesh() {
  if (vtkMesh)
    return vtkMesh;
  else {
    exportToVTKMesh();
    return vtkMesh;
  }
}

void hdf5Reader::exportToMeshBase() {
  std::cout << "Exporting hdf5Reader object to meshBase format..." << std::endl;
  if (!vtkMesh) {
    exportToVTKMesh();
  }
  myMeshBase = meshBase::Create(vtkMesh, outputFname);
}

meshBase *hdf5Reader::getMeshBase() { return myMeshBase; }

void hdf5Reader::setFields(meshBase *myMeshBase,
                           std::vector<std::string> dataNames,
                           std::vector<std::vector<double>> data,
                           int pointOrCell) {
  if (pointOrCell == 0) {
    std::cout << "Setting hdf5Reader object point fields..." << std::endl;
  } else if (pointOrCell == 1) {
    std::cout << "Setting hdf5Reader object cell fields..." << std::endl;
  }

  // Set element data
  auto dataNameItr = dataNames.begin();
  auto dataItr = data.begin();
  while (dataNameItr != dataNames.end() || dataItr != data.end()) {
    std::cout << "    " << *dataNameItr << std::endl;

    if (pointOrCell == 0) {
      myMeshBase->setPointDataArray((*dataNameItr).c_str(), *dataItr);
    } else if (pointOrCell == 1) {
      myMeshBase->setCellDataArray((*dataNameItr).c_str(), *dataItr);
    }

    if (dataNameItr != dataNames.end()) {
      dataNameItr++;
    }
    if (dataItr != data.end()) {
      dataItr++;
    }
  }
}

void hdf5Reader::writeVTK() {
  std::cout << "Writing hdf5Reader object to VTK file" << std::endl;
  if (!myMeshBase) {
    exportToMeshBase();
  }
  myMeshBase->write();
}
