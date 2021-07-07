#include "Transfer/TransferBase.H"

#include <vtkCellData.h>
#include <vtkPointData.h>

 // transfer point data with given names from source to target
 // (converts names to ids before passing to subclass)
int TransferBase::transferPointData(const std::vector<std::string> &arrayNames,
                                    const std::vector<std::string> &newnames)
{
  vtkPointData* pointData = source->getDataSet()->GetPointData();
  std::vector<int> arrayIDs = getArrayIDs(arrayNames, pointData);

  return transferPointData(arrayIDs, newnames);
}

// transfer cell data with given names from source to target - converts names to ids
// (see above)
int TransferBase::transferCellData(const std::vector<std::string> &arrayNames,
                                   const std::vector<std::string> &newnames)
{
  vtkCellData* cellData = source->getDataSet()->GetCellData();
  std::vector<int> arrayIDs = getArrayIDs(arrayNames, cellData);

  return transferCellData(arrayIDs, newnames);
}

std::vector<int> TransferBase::getArrayIDs(const std::vector<std::string>& arrayNames, vtkFieldData* fieldData)
{
  std::vector<int> arrayIDs;
  for(auto arrayName : arrayNames)
  {
    int arrayIndex = getDataArrayIndex(arrayName, fieldData);
    if(arrayIndex == -1)
    {
      std::cerr << "Array " << arrayName << " not found." << std::endl;
      exit(1);
    }
    arrayIDs.push_back(arrayIndex);
  }
  return arrayIDs;
}

int TransferBase::getDataArrayIndex(const std::string& arrayName, vtkFieldData* data)
{
  for(int arrayIndex = 0; arrayIndex < data->GetNumberOfArrays(); ++arrayIndex)
  {
    if(arrayName == data->GetArrayName(arrayIndex)) return arrayIndex;
  }
  return -1;
}

