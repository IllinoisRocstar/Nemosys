#ifndef NEMOSYS_DATASETREGIONSURFACEFILTER_H_
#define NEMOSYS_DATASETREGIONSURFACEFILTER_H_

/*
 * Adopted from vtkDataSetRegionSurfaceFilter and ported from ParaView's
 * VTKExtensions, with slight modifications to deal with extracting curves from
 * 2d meshes as well as maintaining the ID of the twin of elements across an
 * interface.
 */

#include "nemosys_export.h"
#include <vtkCharArray.h>
#include <vtkDataSetSurfaceFilter.h>

namespace NEM {
namespace MSH {

class NEMOSYS_NO_EXPORT dataSetRegionSurfaceFilter
    : public vtkDataSetSurfaceFilter {
 public:
  static dataSetRegionSurfaceFilter *New();
  vtkTypeMacro(dataSetRegionSurfaceFilter, vtkDataSetSurfaceFilter)
  void PrintSelf(ostream &os, vtkIndent indent) override;

  //@{
  /**
   * The name of the cell based array that we use to extract interfaces from
   * Default is "Regions"
   */
  vtkSetStringMacro(RegionArrayName);
  vtkGetStringMacro(RegionArrayName);
  //@}

  int UnstructuredGridExecute(vtkDataSet *input, vtkPolyData *output) override;

  // make it clear we want all the recordOrigCellId signatures from our parent
  using vtkDataSetSurfaceFilter::RecordOrigCellId;

  // override one of the signatures
  void RecordOrigCellId(vtkIdType newIndex, vtkFastGeomQuad *quad) override;

  //@{
  /**
   * Whether to return single sided material interfaces or double sided
   * Default is single
   */
  vtkSetMacro(SingleSided, bool);
  vtkGetMacro(SingleSided, bool);
  //@}

  //@{
  /**
   * The name of the field array that has characteristics of each material
   * Default is "material_properties"
   */
  vtkSetStringMacro(MaterialPropertiesName);
  vtkGetStringMacro(MaterialPropertiesName);
  //@}

  //@{
  /**
   * The name of the field array that has material type identifiers in it
   * Default is "material_ids"
   */
  vtkSetStringMacro(MaterialIDsName);
  vtkGetStringMacro(MaterialIDsName);
  //@}

  //@{
  /**
   * The name of the output field array that records parent materials of each
   * interface Default is "material_ancestors"
   */
  vtkSetStringMacro(MaterialPIDsName);
  vtkGetStringMacro(MaterialPIDsName);
  //@}

  //@{
  /**
   * The name of the field array that has material interface type identifiers in
   * it Default is "interface_ids"
   */
  vtkSetStringMacro(InterfaceIDsName);
  vtkGetStringMacro(InterfaceIDsName);
  //@}

  //@{
  /**
   * The name of the field array that has material interface type identifiers in
   * it Default is "interface_ids"
   */
  vtkSetMacro(Dimension, int);
  vtkGetMacro(Dimension, int);
  //@}

 protected:
  dataSetRegionSurfaceFilter();
  ~dataSetRegionSurfaceFilter() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  /// Implementation of the algorithm.
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  virtual void InsertQuadInHash(vtkIdType a, vtkIdType b, vtkIdType c,
                                vtkIdType d, vtkIdType sourceId,
                                vtkIdType faceId);
  void InsertQuadInHash(vtkIdType a, vtkIdType b, vtkIdType c, vtkIdType d,
                        vtkIdType sourceId) override {
    this->InsertQuadInHash(a, b, c, d, sourceId,
                           -1);  // for -Woverloaded-virtual comp warning
  }

  void InsertTriInHash(vtkIdType a, vtkIdType b, vtkIdType c,
                       vtkIdType sourceId, vtkIdType faceId) override;
  virtual void InsertTriInHash(vtkIdType a, vtkIdType b, vtkIdType c,
                               vtkIdType sourceId) {
    this->InsertTriInHash(a, b, c, sourceId,
                          -1);  // for -Woverloaded-virtual comp warning
  }
  void InsertEdgeInHash(vtkIdType a, vtkIdType b, vtkIdType sourceId,
                        vtkIdType faceId = -1);

  virtual vtkFastGeomQuad *GetNextVisibleQuadFromHash();

  char *RegionArrayName;
  vtkIntArray *RegionArray;
  vtkIdTypeArray *OrigCellIds;
  vtkIntArray *CellFaceIds;
  vtkIdTypeArray *TwinIds;
  bool SingleSided;
  char *MaterialPropertiesName;
  char *MaterialIDsName;
  char *MaterialPIDsName;
  char *InterfaceIDsName;
  int Dimension;

  class Internals;
  Internals *Internal;

 private:
  dataSetRegionSurfaceFilter(const dataSetRegionSurfaceFilter &) = delete;
  void operator=(const dataSetRegionSurfaceFilter &) = delete;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_DATASETREGIONSURFACEFILTER_H_
