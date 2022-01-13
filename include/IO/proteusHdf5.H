/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_PROTEUSHDF5_H_
#define NEMOSYS_PROTEUSHDF5_H_

#include <string>

#include "IO/hdf5Reader.H"
#include "Mesh/meshBase.H"

/* Special purpose class for Proteus HDF5 files */
/** @brief stores information for each block of Proteus data
 **/
struct proteusBlock {
  // block info
  int numElements;           /**< number of elements in block */
  int numVerticesPerElement; /**< number of vertices per element (only one
                               element type allowed per block) */
  std::string blockName;     /**< block name */

  // element type
  int vtkElementType;         /**< VTK element type */
  int originalVtkElementType; /**< original VTK element type. This is different
                                 if we are converting a high-order element into
                                 a low-order element; often useful for
                                 visualization */

  // vertex info
  std::vector<std::vector<std::vector<double>>> vertices; /**< vertices */
  std::vector<int> loc2Glob; /**< map between local and global IDs */

  // field data
  std::vector<std::vector<double>> elementData; /**< element field data */
  std::vector<std::vector<double>> vertexData;  /**< vertex field data */
  proteusBlock(){};
};

/** @brief stores information for each superblock of Proteus data
 **/
struct proteusSuperBlock {
  // global mesh info
  int numElements;              /**< global number of elements */
  int numVertices;              /**< global number of vertices */
  int maxNumVerticesPerElement; /**< maximum number of vertices for all element
                                   types */

  // vertex coordinate info
  std::vector<std::vector<double>>
      coordinates; /**< coordinates[vertex id][dim] = coordinate */
  std::vector<double> xCrd, yCrd,
      zCrd; /**< coordinates separated out into x,y,z coordinates */

  // element info
  std::vector<std::vector<int>>
      elements; /**< connectivity info, elements[id][vertex no] = vertex id */
  std::vector<int>
      elementTypes; // element type info, elementType[id] = vtk element type */

  // element types list
  std::vector<int> elementTypesList; /**< list of element types, ordered by type
                                        according to elements vector */

  std::vector<int>
      originalElementTypesList; /**< list of original element types, ordered by
                                   type according to elements vector */

  // element connectivity sorted by element type and global id
  std::vector<std::map<int, std::vector<int>>>
      vtkConnectivity; /**< connectivity array for all element types by VTK
                          convention, vtkConnectivity[element type][global
                          element id][vertex no] = global vertex id */

  // field data
  std::vector<std::vector<double>>
      elementData; /**< elementData[field no][element id] = data */
  std::vector<std::vector<double>>
      vertexData; /**< cellData[field no][vertex id] = data */

  proteusSuperBlock(){};
};

/** @brief Class to store HDF5 data from Proteus files
 **/
class proteusHdf5 : public hdf5Reader {
public:
  /** @brief Construct Proteus HDF5 object
      @param fieldFName Proteus format HDF5 field filename
      @param meshFName Output VTK mesh filename
      @param edgeSidesetName Name of sideset written to output Exodus file
      @param exoMeshFName Output Exodus mesh filename
      @param lowOrder Boolean converting high order cells to low order
      @param bndryConst Boolean to employ boundary constraint during refinement
   **/
  proteusHdf5(std::string fname, std::string outputFname,
              std::string edgeSidesetName, std::string exoMeshFName,
              bool lowOrder = false, bool bndryConst = true);

  ~proteusHdf5();

  /** @brief Get names of Proteus vector fields
      @param name String description of Proteus vector field names to read
      @param stringVector Vector of Proteus vector field names
      @param numVectors Number of vector fields
   **/
  void getVectorNames(std::string name, std::vector<std::string> &stringVector,
                      int numVectors);

  /** @brief Read information for each Proteus block
      @param group HDF5 Group object containing data
      @param myBlock Proteus block object to read
   **/
  void getBlockInfo(Group &group, proteusBlock &myBlock);

  /** @brief Get coordinates of Proteus block
      @param group HDF5 Group object containing data
      @param myBlock Proteus block object to read
   **/
  void getBlockXYZ(Group &group, proteusBlock &myBlock);

  /** @brief Get global coordinates of Proteus block
      @param group HDF5 Group object containing data
      @param myBlock Proteus block object to read
   **/
  void getBlockGlobalID(Group &group, proteusBlock &myBlock);

  /** @brief Get element data for block
      @param group HDF5 Group object containing data
      @param myBlock Proteus block object to read
   **/
  void getBlockElementData(Group &group, proteusBlock &myBlock);

  /** @brief Get vertex data for block
      @param group Vector of Proteus vector field names
      @param myBlock Number of vector fields
   **/
  void getBlockVertexData(Group &group, proteusBlock &myBlock);

  // access control data
  /**@brief Get number of blocks in Proteus file
     @return Number of blocks
   **/
  int getNumBlocks() const;

  /**@brief Get number of spatial dimensions in Proteus file
     @return Number of dimensions
   **/
  int getNumDimensions() const;

  /**@brief Get number of vertex vectors in Proteus file
     @return Number of vertex vectors
   **/
  int getNumVertexVectors() const;

  /**@brief Get number of element vectors in Proteus file
     @return Number of element vectors
   **/
  int getNumElementVectors() const;

  /**@brief Get character string length of names in Proteus file
     @return Length of string
   **/
  int getCharStringLength() const;

  ///** @brief Write boundary condition information to json file
  //    @param myMeshBase meshBase object of Proteus mesh
  //    @param edgeSidesetName Name of sideset
  // **/
  // void writeBcFile(meshBase *myMeshBase, std::string edgeSidesetName);

  /** @brief Get boundary edge points of mesh using MAdLib
      @param startPts list containing first vertex of each boundary edge
      @param endPts list containing second vertex of each boundary edge
   **/
  void getBoundaryEdgePts(vtkSmartPointer<vtkPoints> startPts,
                          vtkSmartPointer<vtkPoints> endPts);

  /** @brief Get boundary sidesets of mesh
      @param myMeshBase meshBase object of Proteus mesh
      @param startPts list containing first vertex of each boundary edge
      @param endPts list containing second vertex of each boundary edge
      @param sidesetElementIdList list containing element IDs for each edge in
   boundary sideset
      @param sidesetSideIdList list containing side IDs for each edge in
   boundary sideset
   **/
  void getBoundarySideSets(meshBase *myMeshBase,
                           vtkSmartPointer<vtkPoints> startPts,
                           vtkSmartPointer<vtkPoints> endPts,
                           vtkSmartPointer<vtkIdList> sidesetElementIdList,
                           vtkSmartPointer<vtkIdList> sidesetSideIdList);

  /** @brief Get normal vector for 2d meshes
      @param myMeshBase meshBase object of Proteus mesh
      @return normal vector convention, projection onto positive z-axis:
   positive (true) or negative (false)
   **/
  bool get2dCellNodeOrder(meshBase *myMeshBase);

private:
  /** @brief Read Proteus blocks
   **/
  void getBlocks();

  /** @brief Read Proteus control information
   **/
  void getControlInformation();

  /** @brief Merge Proteus block together
   **/
  void mergeBlocks();

  // management data
private:
  bool lowOrder; /**< Boolean converting high order cells to low order (useful
                    for visualization)*/

  // Proteus control data
  int numBlocks;         /**< number of Proteus blocks */
  int numDimensions;     /**< number of spatial dimensions */
  int numVertexVectors;  /**< number of vertex vector fields in Proteus file */
  int numElementVectors; /**< number of element vector fields in Proteus file
                          */
  int charStringLength;  /**< string length of each Proteus variable name */

  // field data names
  std::vector<std::string>
      vertexVectorNames; /**< vector containing names of vertex vectors */
  std::vector<std::string>
      elementVectorNames; /**< vector containing names of element vectors */

  proteusSuperBlock
      mySuperBlock; /**< global block structure for entire Proteus mesh */

  std::vector<proteusBlock>
      proteusBlocks; /**< vector of all individual Proteus-style blocks */
};

#endif  // NEMOSYS_PROTEUSHDF5_H_
