#ifndef TRANSFERBASE_H
#define TRANSFERBASE_H

#include "nemosys_export.h"
#include "meshBase.H"

#include "vtkStaticCellLocator.h"
#include "vtkStaticPointLocator.h"

class NEMOSYS_EXPORT TransferBase {
 public:
  TransferBase()
      : source(nullptr),
        target(nullptr),
        checkQual(false),
        continuous(false),
        c2cTrnsDistTol(1.e-6) {
    std::cout << "TransferBase constructed" << std::endl;
  }

  virtual ~TransferBase() {
    std::cout << "TransferBase destroyed" << std::endl;
  }

 public:
  /**
   * @brief Transfer point data with given ids from source to target.
   * @param arrayIDs array of array ids to specify which fields to transfer
   * @param newnames optional array of names to be applied to transferred fields
   * @return 0 upon completion
   */
  virtual int transferPointData(const std::vector<int> &arrayIDs,
                                const std::vector<std::string> &newnames =
                                    std::vector<std::string>()) = 0;
  /**
   * @brief Transfer cell data with given ids from source to target.
   * @param arrayIDs array of array ids to specify which fields to transfer
   * @param newnames optional array of names to be applied to transferred fields
   * @return 0 upon completion
   */
  virtual int transferCellData(const std::vector<int> &arrayIDs,
                               const std::vector<std::string> &newnames =
                                   std::vector<std::string>()) = 0;
  /**
   * @brief Transfer point data with given field names from source to target.
   * @param arrayIDs array of array names to specify which fields to transfer
   * @param newnames optional array of names to be applied to transferred fields
   * @return 0 upon completion
   */
  int transferPointData(
      const std::vector<std::string> &arrayNames,
      const std::vector<std::string> &newnames = std::vector<std::string>());

  /**
   * @brief Transfer cell data with given field names from source to target.
   * @param arrayIDs array of array names to specify which fields to transfer
   * @param newnames optional array of names to be applied to transferred fields
   * @return 0 upon completion
   */
  int transferCellData(
      const std::vector<std::string> &arrayNames,
      const std::vector<std::string> &newnames = std::vector<std::string>());

  /**
   * @brief Transfer all fields.
   * @param newnames optional array of names to be applied to transferred fields
   * @return 0 upon completion
   */
  virtual int run(const std::vector<std::string> &newnames =
                      std::vector<std::string>()) = 0;

  void setCheckQual(bool x) { checkQual = x; }

  void setContBool(bool x) { continuous = x; }

 protected:
  meshBase *source;
  meshBase *target;

  vtkSmartPointer<vtkStaticCellLocator> srcCellLocator = nullptr;
  vtkSmartPointer<vtkStaticCellLocator> trgCellLocator = nullptr;

  vtkSmartPointer<vtkStaticPointLocator> srcPointLocator = nullptr;
  vtkSmartPointer<vtkStaticPointLocator> trgPointLocator = nullptr;

  bool checkQual;
  bool continuous;  // switch on / off weighted averaging for cell transfer
  double c2cTrnsDistTol;

 private:
  /**
   * @brief given array names and field data, return vector of corresponding
   * array ids in the field data
   */
  std::vector<int> getArrayIDs(const std::vector<std::string> &arrayNames,
                               vtkFieldData *fieldData);
  /**
   * @brief given array name and field data, returns index of array with given
   * name (-1 if not found)
   * @param arrayName name of array
   * @param data mesh field data that includes the sought after arrays
   * @return array id
   */
  int getDataArrayIndex(const std::string &arrayName, vtkFieldData *data);
};

#endif
