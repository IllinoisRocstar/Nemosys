#ifndef NEMOSYS_NUCMESHSHAPEDATA_H_
#define NEMOSYS_NUCMESHSHAPEDATA_H_

#include "nemosys_export.h"
#include "Geometry/ShapeData.H"

#include <string>

#include <SMDSAbs_ElementType.hxx>

#include "Geometry/GeoManager.H"

class SMESH_Gen;
class SMESH_Mesh;
class SMESH_Hypothesis;

namespace NEM {
namespace NUCMESH {

class NEMOSYS_EXPORT NucMeshShapeData : public NEM::GEO::ShapeData {
 public:
  /**
   * Method that is called when a @c GeoManager instance tries to mesh @p shape
   * @details Default implementation does nothing
   * @param[in] shape Shape to mesh
   * @param[in,out] generator Generator that will mesh @shape
   * @param[in,out] mesh Mesh that will mesh @shape
   * @param[out] generatedHyps Output list of created hypotheses/algorithms
   * (for memory management; will be deleted when the @c GeoManager instance is
   * deleted)
   */
  virtual void setupAlgos(
      const TopoDS_Shape &shape, SMESH_Gen &generator, SMESH_Mesh &mesh,
      std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps) const {};
};

/**
 * CRTP pattern to implement @c NEM::GEO::ShapeData::updateModified by copying
 * @tparam T Shape type
 * @tparam Base Immediate parent type
 */
template <typename T, typename Base>
class NEMOSYS_EXPORT CopyOverrideShapeData : public Base {
  using Base::Base;

 public:
  void updateModified(const TopoDS_Shape &shape,
                      const TopTools_ListOfShape &modifiedShapes,
                      NEM::GEO::GeoManager &geoMetadata) const override {
    for (auto &modifiedShape : modifiedShapes) {
      geoMetadata.getMap()[modifiedShape] =
          std::shared_ptr<NEM::GEO::ShapeData>{
              new T{*static_cast<const T *>(this)}};
    }
  }

 protected:
  using CRTPBase = CopyOverrideShapeData<T, Base>;
};

/**
 * Class to hold the material/physical group of a geometry
 */
class NEMOSYS_EXPORT GroupData : public NucMeshShapeData {
 public:
  GroupData(std::string groupName, SMDSAbs_ElementType type)
      : groupName(std::move(groupName)), type(type) {}

  void setupAlgos(const TopoDS_Shape &shape, SMESH_Gen &generator,
                  SMESH_Mesh &mesh,
                  std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps)
      const override;

  std::string groupName;
  SMDSAbs_ElementType type;
};

/**
 *
 */
class NEMOSYS_EXPORT SideSetEdge
    : public CopyOverrideShapeData<SideSetEdge, GroupData> {
 public:
  explicit SideSetEdge(std::string groupName)
      : CRTPBase(std::move(groupName), SMDSAbs_Edge) {}

  void setupAlgos(const TopoDS_Shape &shape, SMESH_Gen &generator,
                  SMESH_Mesh &mesh,
                  std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps)
      const override {
    this->GroupData::setupAlgos(shape, generator, mesh, generatedHyps);
  }
};

/**
 * Class to describe meshing an edge by @c StdMeshers_NumberOfSegments
 */
class NEMOSYS_EXPORT EdgeSegments
    : public CopyOverrideShapeData<EdgeSegments, SideSetEdge> {
 public:
  explicit EdgeSegments(std::string groupName, int numSegments);

  void setupAlgos(const TopoDS_Shape &shape, SMESH_Gen &generator,
                  SMESH_Mesh &mesh,
                  std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps)
      const override;

 private:
  int numSegments_;
};

/**
 * Apply a tri mesh to a face
 */
class NEMOSYS_EXPORT TriMeshSurface
    : public CopyOverrideShapeData<TriMeshSurface, GroupData> {
 public:
  explicit TriMeshSurface(std::string groupName)
      : CRTPBase(std::move(groupName), SMDSAbs_Face) {}

  void setupAlgos(const TopoDS_Shape &shape, SMESH_Gen &generator,
                  SMESH_Mesh &mesh,
                  std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps)
      const override;
};

/**
 * Structured quad mesh (assuming edges meshed appropriately) or quad-dominant
 * mesh
 */
class NEMOSYS_EXPORT QuadMeshSurface
    : public CopyOverrideShapeData<QuadMeshSurface, GroupData> {
 public:
  explicit QuadMeshSurface(std::string groupName)
      : CRTPBase(std::move(groupName), SMDSAbs_Face) {}

  void setupAlgos(const TopoDS_Shape &shape, SMESH_Gen &generator,
                  SMESH_Mesh &mesh,
                  std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps)
      const override;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_NUCMESHSHAPEDATA_H_
