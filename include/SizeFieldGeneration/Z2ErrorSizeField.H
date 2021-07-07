#ifndef NEMOSYS_Z2ERRORSIZEFIELD_H_
#define NEMOSYS_Z2ERRORSIZEFIELD_H_

#include "nemosys_export.h"
#include "SizeFieldGeneration/SizeFieldGen.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT Z2ErrorSizeField : public SizeFieldBase {
  // constructors and destructors
 public:
  Z2ErrorSizeField(vtkDataSet *_ds, int arrayID, int _order);

  ~Z2ErrorSizeField() override {
    std::cout << "Z2ErrorSizeField destroyed" << std::endl;
  }

  // methods
 public:
  // perform patch recovery to get RMSE of recovered values of arrayID
  double computeNodalError(int arrayID) const;
  // compute size field and insert as cell data into mesh's dataSet
  void computeSizeField(vtkDataArray *da) override;

 public:
  int getOrder() const { return order; }
  void setOrder(int _order) { order = _order; }

 private:
  int order;
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_Z2ERRORSIZEFIELD_H_
