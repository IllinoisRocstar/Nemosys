#ifndef VALSIZEFIELD_H
#define VALSIZEFIELD_H

#include "nemosys_export.h"
#include "SizeFieldGen.H"

#include <vector>

class NEMOSYS_EXPORT ValSizeField : public SizeFieldBase
{
  // constructors and destructors
  public:
    ValSizeField(meshBase *_mesh, int arrayID, double _dev_mult,
                 bool _maxIsmin);

    ~ValSizeField() override
    {
      std::cout << "ValSizeField destroyed" << std::endl;
    }

  // methods
  public:
    // computes value of point data at a cell center using average of data
    // at points defining cell
    std::vector<double> computeValAtCell(int cellID, int arrayID) const;
    // compute value of point data at center of each cell
    std::vector<std::vector<double>> computeValAtAllCells(int arrayID) const;
    // compute L2 norm of value of point data at center of each cell
    std::vector<double> computeL2ValAtAllCells(int array) const;
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID) override;

//  private:
//    vtkDataArray *da;
//    std::string sfname;
};

#endif
