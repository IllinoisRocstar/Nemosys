/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
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

