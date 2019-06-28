#ifndef Z2ERRORSIZEFIELD_H
#define Z2ERRORSIZEFIELD_H

#include "nemosys_export.h"
#include "SizeFieldGen.H"

class NEMOSYS_EXPORT Z2ErrorSizeField : public SizeFieldBase
{
  // constructors and destructors
  public:
    Z2ErrorSizeField(meshBase *_mesh, int arrayID);

    ~Z2ErrorSizeField() override
    {
      std::cout << "Z2ErrorSizeField destroyed" << std::endl;
    }

    // methods
  public:
    // perform patch recovery to get RMSE of recovered values of arrayID
    double computeNodalError(int arrayID) const;
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID) override;

  private:
    int order;
};

#endif
