#ifndef VALSIZEFIELD_H
#define VALSIZEFIELD_H

#include "nemosys_export.h"
#include "SizeFieldGen.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT ValSizeField : public SizeFieldBase {
  // constructors and destructors
 public:
  ValSizeField(vtkDataSet *_ds, int arrayID, double _dev_mult, bool _maxIsmin);

  ~ValSizeField() override {
    std::cout << "ValSizeField destroyed" << std::endl;
  }

  // methods
 public:
  // computes value of point data at a cell center using average of data
  // at points defining cell
  static std::vector<double> computeValAtCell(vtkIdList *cell_point_ids,
                                              vtkDataArray *da);
  // compute value of point data at center of each cell
  static std::vector<std::vector<double>> computeValAtAllCells(
      vtkDataSet *ds, vtkDataArray *da);
  // compute L2 norm of value of point data at center of each cell
  static std::vector<double> computeL2ValAtAllCells(vtkDataSet *ds,
                                                    vtkDataArray *da);
  // compute size field and insert as cell data into mesh's dataSet
  void computeSizeField(vtkDataArray *da) override;
};

}  // namespace ADP
}  // namespace NEM

#endif
