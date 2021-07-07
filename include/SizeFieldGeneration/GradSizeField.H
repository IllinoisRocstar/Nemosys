#ifndef NEMOSYS_GRADSIZEFIELD_H_
#define NEMOSYS_GRADSIZEFIELD_H_

#include "nemosys_export.h"
#include "SizeFieldGeneration/SizeFieldGen.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT GradSizeField : public SizeFieldBase {
  // constructors and destructors
 public:
  GradSizeField(vtkDataSet *_ds, int arrayID, double _dev_mult, bool _maxIsmin);

  ~GradSizeField() override {
    std::cout << "GradSizeField destroyed" << std::endl;
  }

  // methods
 public:
  // computes the gradient of point data at a cell using
  // derivatives of shape interpolation functions
  static std::vector<double> computeGradAtCell(vtkCell *cell, vtkDataArray *da);
  // compute L2 norm of gradient of point data at each cell
  static std::vector<double> computeL2GradAtAllCells(vtkDataSet *ds,
                                                     vtkDataArray *da);
  // compute size field and insert as cell data into mesh's dataSet
  void computeSizeField(vtkDataArray *da) override;
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_GRADSIZEFIELD_H_
