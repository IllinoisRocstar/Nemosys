#ifndef MESHGEN_H
#define MESHGEN_H

#include "nemosys_export.h"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>

#include <string>
#include <iostream>


class meshingParams;

class NEMOSYS_EXPORT meshGen {
  public:
    meshGen() : dataSet(nullptr)
    { std::cout << "Base mesh generator created" << std::endl; }

    virtual ~meshGen()
    { std::cout << "Base mesh generator destroyed" << std::endl; }

    // creates generator with default parameters
    static meshGen *Create(const std::string &fname,
                           const std::string &meshEngine);
    // creates generater with specified parameters
    static meshGen *Create(const std::string &fname,
                           const std::string &meshEngine,
                           meshingParams *params);
    virtual int createMeshFromSTL(const char *fname) = 0;
    vtkSmartPointer<vtkDataSet> getDataSet() const { return dataSet; }

  protected:
    vtkSmartPointer<vtkDataSet> dataSet;
};

#endif
