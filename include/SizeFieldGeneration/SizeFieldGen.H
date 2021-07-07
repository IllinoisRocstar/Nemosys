#ifndef NEMOSYS_SIZEFIELDGEN_H_
#define NEMOSYS_SIZEFIELDGEN_H_

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

#include "nemosys_export.h"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT SizeFieldBase {
  // constructors and destructors
 public:
  SizeFieldBase() : ds(nullptr), dev_mult(1.5), maxIsmin(true), sizeFactor(1.) {
    std::cout << "SizeFieldBase constructed" << std::endl;
  }
  SizeFieldBase(vtkDataSet *_ds, int arrayID, double _dev_mult, bool _maxIsmin,
                const std::string &arrName);

  virtual ~SizeFieldBase() {
    std::cout << "SizeFieldBase destroyed" << std::endl;
  }

  // SizeFieldBase factory methods
 public:
  // will return derived type GradSizeField|ValSizeField|ErrorSF with SF in
  // mesh's dataSet
  //    static SizeFieldBase *
  //    Create(meshBase *_mesh, const std::string &method, int arrayID);
  static SizeFieldBase *Create(vtkDataSet *_dataSet, const std::string &method,
                               int arrayID, double _dev_mult, bool _maxIsmin,
                               double _sizeFactor = 1.0, int _order = 1);
  static std::unique_ptr<SizeFieldBase> CreateUnique(
      vtkDataSet *_dataSet, const std::string &method, int arrayID,
      double _dev_mult, bool _maxIsmin, double _sizeFactor = 1.0,
      int _order = 1);

  // method
 public:
  // compute size field and insert as cell data into mesh's dataSet
  virtual void computeSizeField(vtkDataArray *da) = 0;

  void setSizeFactor(double sf) { sizeFactor = sf; }

 protected:
  vtkSmartPointer<vtkDataSet> ds;
  double dev_mult;
  bool maxIsmin;
  vtkSmartPointer<vtkDataArray> da;
  std::string sfname;
  double sizeFactor;

  // helpers
 protected:
  // identifies cells to refine and mutates current size values
  // into a compatible size field for the mesh
  void mutateValues(std::vector<double> &values) const;
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_SIZEFIELDGEN_H_
