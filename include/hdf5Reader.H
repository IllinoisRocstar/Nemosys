#ifndef __HDF5READER_H
#define __HDF5READER_H

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

// vtk
#include <vtkCellTypes.h>
#include <vtkDataSet.h>
#include <vtkIdList.h>
#include <vtkPoints2D.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// HDF5
#include <H5Cpp.h>
#include <H5Exception.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

class meshBase;

/** @brief Class for reading data from HDF5 files
 **/
class hdf5Reader {
 public:
  /** @brief Construct HDF5 reader object
      @param fname HDF5 filename
      @param _outputFName Filename for VTK output file
      @param _verb Boolean verbosity flag
   **/
  hdf5Reader(std::string fname, std::string _outputFname, int _verb = 0)
      :  // Initialization of class variables
        h5FileName(fname),
        outputFname(_outputFname),
        verb(_verb),
        vtkMesh(0) {
    // Parse filename
    std::size_t _loc = h5FileName.find_last_of(".");
    baseH5FileName = h5FileName.substr(0, _loc);
  };

  // Destructor
  virtual ~hdf5Reader(){};

  // File I/O
  /** @brief Open HDF5 file
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  int openFile();

  // File I/O
  /** @brief Close HDF5 file
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  int closeFile();

  // Read data
  //    Reading string data
  // File I/O
  /** @brief Read string data from top level HDF5 DataSet
      @param dsetName Dataset name
      @param buffer Buffer into which string data is read
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  int readTopDataSet(std::string dsetName, std::string &buffer);
  /** @brief Read string data from generic HDF5 DataSet
      @param dset HDF5 dataset
      @param buffer Buffer into which data is read
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  int readDataSet(DataSet dset, std::string &buffer);
  //    Reading numeric data
  /** @brief Read numeric data from top level HDF5 DataSet
      @param dsetName Dataset name
      @param buffer Buffer into which numeric data is read
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  template <class T>
  int readTopDataSet(std::string dsetName, std::vector<T> &buffer) {
    try {
      std::cout << "H5_VERS_MAJOR = " << H5_VERS_MAJOR << std::endl;
      std::cout << "H5_VERS_MINOR = " << H5_VERS_MINOR << std::endl;
      std::cout << "H5_VERS_RELEASE = " << H5_VERS_RELEASE << std::endl;
      std::cout << "H5_VERS_SUBRELEASE = " << H5_VERS_SUBRELEASE << std::endl;
      // Open dataset
      DataSet dset = file.openDataSet(dsetName);
      readDataSet(dset, buffer);
    }
    // catch failure caused by the H5File operations
    catch (FileIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
      error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }

    return 0;
  }

  //    Reading numeric data
  /** @brief Read numeric data from generic HDF5 DataSet
      @param dset HDF5 dataset
      @param buffer Buffer into which numeric data is read
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  template <class T> int readDataSet(DataSet dset, std::vector<T> &buffer) {
    try {
      // Get type of data
      DataType dtype = dset.getDataType();
      // Read data
      dset.read(&buffer[0], dtype);
      // Close dataset
      dset.close();
    }
    // catch failure caused by the H5File operations
    catch (DataSetIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
      error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }

    return 0;
  }

  // Read groups
  /** @brief Read existing HDF5 Group
      @param groupName HDF5 group object name
      @param group HDF5 group object
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  int readGroup(std::string groupName, Group &group);

  // Reading numeric data from group
  /** @brief Read numeric data from HDF5 group
      @param dsetName DataSet name
      @param group HDF5 group object
      @param buffer Buffer into which numeric data is read
      @return Integer flag: 0 (failure) or 1 (success)
   **/
  template <class T>
  int readGroupDataSet(std::string dsetName, Group &group,
                       std::vector<T> &buffer) {
    try {
      // Open dataset
      DataSet dset = group.openDataSet(dsetName);
      // Get type of data
      DataType dtype = dset.getDataType();
      // Read data
      dset.read(&buffer[0], dtype);
      // Close dataset
      dset.close();
    }
    // catch failure caused by the H5File operations
    catch (GroupIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
      error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    } catch (DataSetIException error) {
#if H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 8 && H5_VERS_RELEASE >= 20
      error.printErrorStack();
#elif H5_VERS_MAJOR >= 1 && H5_VERS_MINOR >= 10 && H5_VERS_RELEASE >= 2
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }

    return 0;
  }

  // Getter functions
  /** @brief Get HDF5 filename
    @return HDF5 filename
 **/
  std::string getFileName();

  /** @brief Flatten 2d data buffer
      @param originalBuffer Original unflattened data
      @param flatBuffer Output flattened data
  **/
  template <class T>
  void flattenBuffer2d(std::vector<std::vector<T>> originalBuffer,
                       std::vector<T> &flatBuffer) {
    // Populate flat buffer
    for (auto bufItr1 = originalBuffer.begin(); bufItr1 != originalBuffer.end();
         bufItr1++) {
      for (auto bufItr2 = (*bufItr1).begin(); bufItr2 != (*bufItr1).end();
           bufItr2++) {
        flatBuffer.push_back(*bufItr2);
      }
    }
  }

  /** @brief Unflatten 2d data buffer
      @param flatBuffer Original flattened data
      @param originalBuffer Output unflattened data
  **/
  template <class T>
  void unflattenBuffer2d(std::vector<T> flatBuffer,
                         std::vector<std::vector<T>> &originalBuffer) {
    // Populate original buffer
    auto flatItr = flatBuffer.begin();
    for (auto bufItr1 = originalBuffer.begin(); bufItr1 != originalBuffer.end();
         bufItr1++) {
      for (auto bufItr2 = (*bufItr1).begin(); bufItr2 != (*bufItr1).end();
           bufItr2++) {
        *bufItr2 = *flatItr;
        flatItr++;
      }
    }
  }

  // General implementation for flattening of arbitrary type and dimension
  // vector
  // template <template<typename...> class R=std::vector,
  //          typename Top,
  //          typename Sub = typename Top::value_type>
  //   R<typename Sub::value_type> flattenBuffer(Top const& all)
  //{
  //    using std::begin;
  //    using std::end;
  //
  //    R<typename Sub::value_type> accum;
  //
  //    for(auto& sub : all)
  //        std::copy(begin(sub), end(sub), std::inserter(accum, end(accum)));
  //
  //    return accum;
  //}

  // export to vtk format without going through madlib/gmsh (needed for
  // simmetrix interface)
  /** @brief Export to VTK format
   **/
  void exportToVTKMesh();
  /** @brief Get VTK mesh object
      @return Pointer to vtkDataSet
  **/
  vtkSmartPointer<vtkDataSet> getVTKMesh();
  /** @brief Export HDF5 mesh to meshBase object
   **/
  void exportToMeshBase();
  /** @brief Get meshBase object
      @return meshBase object
  **/
  meshBase *getMeshBase();
  /** @brief Write mesh to VTK file
   **/
  void writeVTK();

  // Setting meshBase fields
  /** @brief Set meshBase field data
      @param myMeshBase meshBase object
      @param dataNames Vector of names of data fields
      @param data Vector storing data
      @param pointOrCell Boolean indicating point or cell data
  **/
  void setFields(meshBase *myMeshBase, std::vector<std::string> dataNames,
                 std::vector<std::vector<double>> data, int pointOrCell);

  // Setter functions
  // Setting meshBase fields
  /** @brief Set number of vertices in meshBase
      @param numVertices Number of vertices
  **/
  void setNumberOfVertices(int numVertices);
  /** @brief Set number of elements in meshBase
      @param numElements Number of elements
  **/
  void setNumberOfElements(int numElements);
  /** @brief Set number of spatial dimensions in meshBase
      @param numDimensions Number of spatial dimensions
  **/
  void setNumberOfDimensions(int numDimensions);
  /** @brief Set meshBase 1d coordinate data
      @param xCrd vector of x-coordinates
  **/
  void setCoordinates(std::vector<double> xCrd);
  /** @brief Set meshBase 2d coordinate data
      @param xCrd vector of x-coordinates
      @param yCrd vector of y-coordinates
  **/
  void setCoordinates(std::vector<double> xCrd, std::vector<double> yCrd);
  /** @brief Set meshBase 3d coordinate data
      @param xCrd vector of x-coordinates
      @param yCrd vector of y-coordinates
      @param zCrd vector of z-coordinates
  **/
  void setCoordinates(std::vector<double> xCrd, std::vector<double> yCrd,
                      std::vector<double> zCrd);
  /** @brief Set meshBase element connectivity table
      @param elementConnectivity Element connectivity table
  **/
  void setConnectivities(
      std::vector<std::map<int, std::vector<int>>> elementConnectivity);
  /** @brief Set VTK element types for each element
      @param vtkElementTypes Vector of VTK element types for each element
  **/
  void setElementTypes(std::vector<int> vtkElementTypes);
  /** @brief Set list of each unique type of VTK element in mesh
      @param vtkElementTypesList Vector of each unique VTK element type in mesh
  **/
  void setElementTypesList(std::vector<int> vtkElementTypesList);

  /** @brief Get element connectivity for global element ID
      @param elementId Element ID
      @return Element connectivity table
  **/
  std::vector<int> getElementConnectivity(int elementId);

 protected:
  // class variables
  std::string h5FileName;              /**< HDF5 filename */
  std::string baseH5FileName;          /**< base of HDF5 filename */
  std::string outputFname;             /**< VTK output filename */
  int verb;                            /**< verbosity flag */
  H5File file;                         /**< HDF5 file object */
  vtkSmartPointer<vtkDataSet> vtkMesh; /**< VTK DataSet object */
  meshBase *myMeshBase;                /**< meshBase object */

  // mesh data
  std::vector<double> xCrd, yCrd, zCrd; /**< vector of coordinates */
  std::vector<std::map<int, std::vector<int>>>
      elementConnectivity; /**< element connectivity table */
  int numVertices;         /**< number of vertices in mesh */
  int numElements;         /**< number of elements in mesh */
  int numDimensions;       /**< number of spatial dimensions */
  std::vector<int>
      vtkElementTypes; /**< vector of VTK element type for each element */
  std::vector<int> vtkElementTypesList; /**< vector storing list of each unique
                                           VTK element type */
};

#endif  // __HDF5READER_H
