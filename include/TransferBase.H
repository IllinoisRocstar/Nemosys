#ifndef TRANSFERBASE_H
#define TRANSFERBASE_H

#include "nemosys_export.h"
#include "meshBase.H"

#include "vtkStaticPointLocator.h"
#include "vtkStaticCellLocator.h"

class NEMOSYS_EXPORT TransferBase
{
  // constructors and destructors
  public:
    TransferBase()
        : source(nullptr), target(nullptr),
          checkQual(false), continuous(false), c2cTrnsDistTol(1.e-6)
    {
      std::cout << "TransferBase constructed" << std::endl;
    }

    virtual ~TransferBase()
    {
      std::cout << "TransferBase destroyed" << std::endl;
    }

  // transfer methods
  public:

    // transfer point data with given id from source to target
    virtual int transferPointData(const std::vector<int> &arrayIDs,
                                  const std::vector<std::string> &newnames = std::vector<std::string>()) = 0;
    // transfer cell data with given id from source to target
    virtual int transferCellData(const std::vector<int> &arrayIDs,
                                 const std::vector<std::string> &newnames = std::vector<std::string>()) = 0;

     // transfer point data with given names from source to target
     // (converts names to ids before passing to subclass)
    int transferPointData(const std::vector<std::string> &arrayNames,
                          const std::vector<std::string> &newnames = std::vector<std::string>());

    // transfer cell data with given names from source to target - converts names to ids
    // (see above)
    int transferCellData(const std::vector<std::string> &arrayNames,
                         const std::vector<std::string> &newnames = std::vector<std::string>());

    // transfer all data from source to target
    virtual int
    run(const std::vector<std::string> &newnames = std::vector<std::string>()) = 0;

    // set whether to check transfer quality
    void setCheckQual(bool x) { checkQual = x; }

    void setContBool(bool x) { continuous = x; }

  protected:
    meshBase *source;
    meshBase *target;

    // default locators, still used by transferCellData
    vtkSmartPointer<vtkStaticCellLocator> srcCellLocator = nullptr;
    vtkSmartPointer<vtkStaticCellLocator> trgCellLocator = nullptr;

    // thread safe locators
    vtkSmartPointer<vtkStaticPointLocator> srcPointLocator = nullptr;
    vtkSmartPointer<vtkStaticPointLocator> trgPointLocator = nullptr;

    bool checkQual;
    bool continuous; // switch on / off weighted averaging for cell transfer
    double c2cTrnsDistTol;

  private:

    // given array names and field data, return vector of corresponding array ids in the field data
    std::vector<int> getArrayIDs(const std::vector<std::string>& arrayNames, vtkFieldData* fieldData);
    // given array name and field data, returns index of array with given name (-1 if not found)
    int getDataArrayIndex(const std::string& arrayName, vtkFieldData* data);


};

#endif
