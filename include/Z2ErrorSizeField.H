#ifndef Z2ERRORSIZEFIELD_H
#define Z2ERRORSIZEFIELD_H

#include <SizeFieldBase.H>

class Z2ErrorSizeField: public SizeFieldBase
{
  // constructors and destructors
  public:
    Z2ErrorSizeField(meshBase* _mesh, int arrayID, int order, double _dev_mult, bool _maxIsmin);
    ~Z2ErrorSizeField()
    {
      std::cout << "Z2ErrorSizeField destroyed" << std::endl;
    }
    
  // methods
  public:
    // perform patch recovery to get rmse of recovered values of arrayID
    double computeNodalError(int arrayID); 
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID);

  private:
    int order;
};



#endif
