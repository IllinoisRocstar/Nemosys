#ifndef VALSIZEFIELD_H
#define VALSIZEFIELD_H

#include <SizeFieldBase.H>

class ValSizeField: public SizeFieldBase
{

  // constructors and destructors
  public:
    ValSizeField(meshBase* _mesh, int arrayID, double _dev_mult, bool _maxIsmin);
    ~ValSizeField()
    {
      std::cout << "ValSizeField destroyed" << std::endl;
    }

  // methods
  public:
    // computes value of point data at a cell center using average of data
    // at points defining cell
    std::vector<double> computeValAtCell(int cellID, int arrayID);
    // compute value of point data at center of each cell
    std::vector<std::vector<double>> computeValAtAllCells(int arrayID);
    // compute L2 norm of value of point data at center of each cell
    std::vector<double> computeL2ValAtAllCells(int array);
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID);

//  private:
//    vtkDataArray* da;
//    std::string sfname;
};

#endif
