#ifndef GRADSIZEFIELD_H
#define GRADSIZEFIELD_H

#include <SizeFieldBase.H>

class GradSizeField: public SizeFieldBase
{
  // constructors and destructors
  public:
    GradSizeField(meshBase* _mesh, int arrayID, double _dev_mult, bool _maxIsmin);
    ~GradSizeField()
    {
      std::cout << "GradSizeField destroyed" << std::endl;
    }
    
  // methods
  public:
    // computes the gradient of point data at a cell using 
    // derivatives of shape interpolation functions
    std::vector<double> computeGradAtCell(int cellID, int arrayID);
    // compute L2 norm of gradient of point data at each cell 
    std::vector<double> computeL2GradAtAllCells(int arrayID);
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID);

//  private:
//    vtkDataArray* da;
//    std::string sfname; 
};



#endif
