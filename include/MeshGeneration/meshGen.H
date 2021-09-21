#ifndef NEMOSYS_MESHGEN_H_
#define NEMOSYS_MESHGEN_H_

#include "nemosys_export.h"
#include "Mesh/geoMeshFactory.H"
#include "Mesh/geoMeshBase.H"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>

#include <string>
#include <iostream>

class NEMOSYS_EXPORT meshGen {
 public:
  meshGen() : dataSet(nullptr) {
    std::cout << "Base mesh generator created" << std::endl;
  }

  virtual ~meshGen() {
    std::cout << "Base mesh generator destroyed" << std::endl;
  }

  virtual int createMeshFromSTL(const char *fname) = 0;
  vtkSmartPointer<vtkDataSet> getDataSet() const { return dataSet; }
  std::unique_ptr<NEM::MSH::geoMeshBase> gmData;

 protected:
  vtkSmartPointer<vtkDataSet> dataSet;
};

#endif  // NEMOSYS_MESHGEN_H_
