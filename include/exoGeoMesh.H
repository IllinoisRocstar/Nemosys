#ifndef NEMOSYS_EXOGEOMESH_H_
#define NEMOSYS_EXOGEOMESH_H_

#include <vtkExodusIIReader.h>
#include <vtkMultiBlockDataSet.h>
#include <map>
#include <set>
#include <vector>

#include "geoMeshBase.H"

namespace NEM {
namespace MSH {

/**
 * @class exoGeoMesh
 * @brief A concrete implementation of @c geoMeshBase representing an Exodus II
 * mesh
 * @details Note that only the title, dimension, element blocks, element block
 * properties, element block variables, node variables, global variables, node
 * sets, and side sets (and the names of each of these fields) are read and
 * stored in an exoGeoMesh class. There is one possible exception: to
 * represent a large number of side sets, if there is a side set variable
 * named @c "NEM_SIDE_SET_ID" then sides are assigned to separate side sets
 * with side set IDs set according to this variable (and all side set names
 * are ignored). Side set sides are ignored if the source element does not have
 * the same dimensionality as the mesh.
 */
class NEMOSYS_EXPORT exoGeoMesh : public geoMeshBase {
 public:
  static exoGeoMesh *New();
  vtkTypeMacro(exoGeoMesh, geoMeshBase)

 public:
  /**
   * Reads an EXODUS II file into an exoGeoMesh object.
   * @param fileName Name of file to read. Empty string interpreted as creating
   * new exoGeoMesh with empty mesh
   * @param timeStep Time step (using EXODUS II indexing) to read
   * @return Pointer to a new exoGeoMesh object
   */
  static exoGeoMesh *Read(const std::string &fileName, int timeStep = 1);

 public:
  exoGeoMesh();

  /**
   * @brief exoGeoMesh constructor from file.
   * @param fileName Name of EXODUS II file to read
   * @param timeStep Time step (using EXODUS II indexing) to read
   */
  explicit exoGeoMesh(const std::string &fileName, int timeStep = 1);

  ~exoGeoMesh() override;

 protected:
  /**
   * @brief exoGeoMesh constructor from vtkExodusIIReader
   * @details Note that side set side indices in the reader's output are ignored
   * and are read from @c reader->GetFileName() instead. Element block
   * properties are also read from the file separately.
   * @param reader Existing vtkExodusIIReader from which to create exoGeoMesh
   * object. Update() should already have been called.
   */
  explicit exoGeoMesh(const vtkSmartPointer<vtkExodusIIReader> &reader);

 public:
  /**
   * @brief Write exoGeoMesh out to a file.
   * @details If @p fileName already exists, it is overwritten. If there are
   * more than 900 side sets, they are concatenated, with a side set variable
   * named @c "NEM_SIDE_SET_ID" added that contains the id of the side in this
   * exoGeoMesh.
   * @param fileName Name of file. Must end in ".exo", ".e", ".gen", or ".g"
   */
  void write(const std::string &fileName) override;
  /**
   * Write a summary of this exoGeoMesh to @p out.
   * @param out Stream to write summary to.
   */
  void report(std::ostream &out) const override;
  /**
   * @copydoc geoMeshBase::reconstructGeo()
   * If @c physGrpName is not empty, then physical groups are set based on the
   * property @c physGrpName for each element block. Otherwise, each element
   * block represents its own physical group.
   */
  void reconstructGeo() override;
  void takeGeoMesh(geoMeshBase *otherGeoMesh) override;

  /**
   * Get the title of the exoGeoMesh (used for writing out to EXODUS II file)
   * @return reference to title
   */
  const std::string &getTitle() const;
  /**
   * Set the title of the exoGeoMesh (used for writing out to EXODUS II file)
   * @param title new title to set
   */
  void setTitle(const std::string &title);
  /**
   * Get number of element blocks
   * @return number of element blocks
   */
  int getNumElemBlocks() const;
  /**
   * Get all element block names
   * @return element block names, indexed by ID
   */
  std::map<int, std::string> getElemBlockNames() const;
  /**
   * Get the name of an element block (used for writing out to EXODUS II file)
   * @param id ID of the element block
   * @return reference to block name
   */
  const std::string &getElemBlockName(int id) const;
  /**
   * Set the name of an element block (used for writing out to EXODUS II file)
   * @param id ID of the element block
   * @param name new name of the element block
   */
  void setElemBlockName(int id, const std::string &name);
  /**
   * Get all element block IDs.
   * @return All element block IDs.
   */
  std::vector<int> getElemBlockIds() const;
  /**
   * Get the element block ID of a cell
   * @param cellIdx 0-based index of a cell
   * @return block ID that the cell belongs to; -1 if not found
   */
  int getElemBlockId(vtkIdType cellIdx) const;
  /**
   * Get (0-based) indices of the cells that are in a specific element block
   * @param blockId ID of an element block
   * @return Indices of all cells that belong to @p blockID
   */
  std::vector<vtkIdType> getElemBlock(int blockId) const;
  /**
   * Add a new element block to the mesh. Invalidates geometry, if any.
   * PointData and CellData arrays are intersected.
   * @param elements vtkUnstructuredGrid of cells to add to the mesh. All
   * cells should be of same cell type.
   * @param name Name of new element block
   * @param tol Tolerance for merging points
   * @return The ID of the new element block.
   */
  int addElemBlock(vtkUnstructuredGrid *elements, const std::string &name = "",
                   float tol = 1e-15f);
  /**
   * Add cells to an existing element block. Invalidates geometry, if any.
   * PointData and CellData arrays are intersected.
   * @param elements vtkUnstructuredGrid of cells to add to the mesh. All
   * cells should be of the same type as cells in the block given by @p blockId.
   * @param blockId ID of existing element block
   * @param tol Tolerance for merging points
   */
  void addCellsToBlock(vtkUnstructuredGrid *elements, int blockId,
                       float tol = 1e-15f);
  /**
   * @brief Reassign cells to a different element block
   * @details Cells must all have same cell type. May result in empty blocks.
   * Invalidates geometry, if any.
   * @param cells vector of cell indices
   * @param blockId id of block to assign to; if not an existing block, one is
   * created (if @p blockId is positive, the new block will have that ID;
   * otherwise, one is chosen)
   * @param elemBlockName if blockId is not an existing block, the name of the
   * new block
   * @return element block that the cells were assigned to
   */
  int reassignCells(const std::vector<vtkIdType> &cells, int blockId = -1,
                    const std::string &elemBlockName = {});
  /**
   * Add a property to element blocks (with default value 0 on all element
   * blocks)
   * @param propName Name of new property
   * @return Whether or not new property was added
   */
  bool addElemBlockProperty(const std::string &propName);
  /**
   * Get the names of all block properties
   * @return Names of all element block properties
   */
  std::vector<std::string> getElemBlockPropertyNames() const;
  /**
   * Get the value of an element block property on a block
   * @param propName Name of an element block property
   * @param blockId ID of an element block
   * @return Value of a property on an element block
   */
  int getElemBlockProperty(const std::string &propName, int blockId) const;
  /**
   * Set the value of an element block property on a block
   * @param propName Name of an element block property
   * @param blockId ID of an element block
   * @param value Value of a property on an element block
   */
  void setElemBlockProperty(const std::string &propName, int blockId,
                            int value);
  /**
   * Get number of side sets
   * @return number of side sets
   */
  int getNumSideSets() const;
  /**
  * Get all side set names
  * @return side set names, indexed by ID
  */
  const std::map<int, std::string> &getSideSetNames() const;
  /**
   * Get the name of a side set (used for writing out to EXODUS II file)
   * @param id ID of a side set
   * @return Reference to the name of a side set
   */
  const std::string &getSideSetName(int id) const;
  /**
   * Set the name of a side set (used for writing out to EXODUS II file)
   * @param id ID of a side set
   * @param name New name of the side set
   */
  void setSideSetName(int id, const std::string &name);
  /**
   * Get all side set IDs.
   * @return All side set IDs.
   */
  std::vector<int> getSideSetIds() const;
  /**
   * Get sides in a side set
   * @param sideSetId ID of a side set
   * @return Pair of vectors, first of which contains the (0-based) index
   * of the cell the side belongs to and the second contains the side (according
   * to VTK face/edge ordering)
   */
  std::pair<std::vector<vtkIdType>, std::vector<int>> getSideSet(
      int sideSetId) const;
  /**
   * Add a new side set
   * @param elements The (0-based) index of the cells that each side comes from
   * @param sides The index of the face/edge that of the side (according to VTK
   * face/edge ordering)
   * @param name Name of new side set
   * @return ID of new side set
   */
  int addSideSet(const std::vector<vtkIdType> &elements,
                 const std::vector<int> &sides, const std::string &name = "");
  /**
   * Get number of node sets
   * @return number of node sets
   */
  int getNumNodeSets() const;
  /**
  * Get all node set names
  * @return node set names, indexed by ID
  */
  std::map<int, std::string> getNodeSetNames() const;
  /**
   * Get the name of a node set (used for writing out to EXODUS II file)
   * @param id ID of a node set
   * @return Reference to the name of a node set
   */
  const std::string &getNodeSetName(int id) const;
  /**
   * Set the name of a node set (used for writing out to EXODUS II file)
   * @param id ID of a node set
   * @param name New name of the node set
   */
  void setNodeSetName(int id, const std::string &name);
  /**
   * Get the IDs of all node sets
   * @return All node set IDs
   */
  std::vector<int> getNodeSetIds() const;
  /**
   * Get the nodes of a node set.
   * @param nodeSetId ID of a node set
   * @return Vector of (0-based) node indices of the nodes in the node set.
   */
  const std::vector<vtkIdType> &getNodeSet(int nodeSetId) const;
  /**
   * Add a new node set
   * @param nodes Node indices that should belong on the new node set (0-based
   * indexing)
   * @param name Name of new node set
   * @return ID of new node set
   */
  int addNodeSet(const std::vector<vtkIdType> &nodes,
                 const std::string &name = "");
  /**
   * Get the name of the element block property that maps element blocks to
   * physical groups (if empty, each element block represents its own physical
   * group)
   * @return Name of element block property mapping element blocks to physical
   * groups (or empty string if not set)
   */
  const std::string &getPhysGrpPropertyName() const;
  /**
   * Sets the name of the element block property that maps element blocks to
   * physical groups (if empty, each element block represents its own physical
   * group)
   * @param physGrpName Name of element block property mapping element blocks to
   * physical groups (use empty string to unset). Must be an existing element
   * block property name.
   */
  void setPhysGrpPropertyName(const std::string &physGrpName);
  /**
   * Stitch another exoGeoMesh object onto this one. PointData and CellData
   * arrays are intersected. Adds side sets and node sets if available. Element
   * block properties from @p otherGM are copied if property also defined in
   * this exoGeoMesh.
   * @param otherGM Other exoGeoMesh to addd to this one
   * @param tol Tolerance for merging points
   */
  void stitch(const exoGeoMesh &otherGM, float tol = 1e-15f);
  /**
   * Scale point locations. Does not change point data.
   * @param scale Multiplicative factor for point locations
   */
  void scaleNodes(double scale);

 private:
  static vtkSmartPointer<vtkExodusIIReader> getExoReader(
      const std::string &fileName, int timeStep = 1);
  static GeoMesh exoReader2GM(vtkSmartPointer<vtkExodusIIReader> reader);
  static const char *getExoIdArrName();
  constexpr static const char *SIDE_SET_ID_VAR_NAME = "NEM_SIDE_SET_ID";
  /**
   * Clears element block names, element block properties, and dimension.
   * Creates element blocks and sets the dimension based on the @c GeoMesh.
   */
  void resetElemBlocksAndDim();
  /**
   * Renumbers the nodes in all node sets based on @p nodeMap if provided.
   * Removes duplicate node ids within a node set.
   * @param nodeMap map from old to new node indices (can be many-to-one)
   */
  void resetNodeSetPoints(vtkIdTypeArray *nodeMap = nullptr);
  /**
   * Renumbers the side set points based on @p nodeMap.
   * @param nodeMap map from old to new node indices (can be many-to-one)
   */
  void resetSideSetPoints(vtkIdTypeArray *nodeMap);
  /**
   * Sets the (Exodus) side set IDs of each cell in the @c GeoMesh.sideSet by
   * creating a side set for each pair (side geometric entity, parent
   * physical group).
   */
  void setSideSetObjId();

  void resetNative() override;

 protected:
  struct elemBlock {
    std::string name;
    VTKCellType cellType;
    std::map<std::string, int> properties;  // missing values treated as 0
  };
  struct nodeSet {
    std::string name;
    std::vector<vtkIdType> nodes;  // correspond to indices in the _geoMesh.mesh
  };

 private:
  std::string _title;
  std::map<int, elemBlock> _elemBlocks;
  std::map<int, std::string> _sideSetNames;
  std::map<int, nodeSet> _nodeSets;
  std::set<std::string> _elemBlockPropNames;
  std::string _physGrpName;  // Should correspond to a property
  int _dim;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_EXOGEOMESH_H_
