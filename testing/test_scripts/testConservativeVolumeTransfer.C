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
#include <Mesh/meshBase.H>
#include <Mesh/vtkMesh.H>
#include <gtest/gtest.h>
#include <Drivers/TransferDriver.H>
#include <Transfer/ConservativeVolumeTransfer.H>
#include <Transfer/ConservativeSurfaceTransfer.H>
#include <Mesh/vtkMesh.H>

#include <chrono>

#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkXMLUnstructuredGridWriter.h"

const char* sourceFile;
const char* targetFile;
const char* outputFile;

vtkMesh * sourceMesh;
vtkMesh * targetMesh;
ConservativeVolumeTransfer * volumeTransfer;

std::chrono::time_point start_time;

TEST(ConservativeVolumeTransferTest, driverTest)
{
  TransferDriver* transferDriver = new TransferDriver(std::string(sourceFile),
                                                      std::string(targetFile),
                                                      std::string("Conservative Volume Transfer"),
                                                      std::string(outputFile),
                                                      false);
  delete transferDriver;
}

TEST(ConservativeVolumeTransferTest, constructionTest)
{
  start_time = std::chrono::system_clock::now();
  sourceMesh = new vtkMesh(sourceFile);
  targetMesh = new vtkMesh(targetFile);

  std::cout << "source mesh has :" << std::endl;
  std::cout << sourceMesh->getDataSet()->GetNumberOfPoints() << " nodes" << std::endl;
  std::cout << sourceMesh->getDataSet()->GetNumberOfCells() << " cells" << std::endl;
  std::cout << std::endl;
  std::cout << "target mesh has " << std::endl;
  std::cout << targetMesh->getDataSet()->GetNumberOfPoints() << " nodes" << std::endl;
  std::cout << targetMesh->getDataSet()->GetNumberOfCells() << " cells" << std::endl;
  std::cout << std::endl;

  volumeTransfer = new ConservativeVolumeTransfer(sourceMesh, targetMesh);
  volumeTransfer->enableSurfaceTransfer();
}

TEST(ConservativeVolumeTransferTest, supermeshConstructionTest)
{
  volumeTransfer->constructSupermesh();
}

TEST(ConservativeVolumeTransferTest, massMatrixTest)
{
  volumeTransfer->constructMassMatrix();
}

TEST(ConservativeVolumeTransferTest, mixedMassMatrixTest)
{
  volumeTransfer->constructMixedMassMatrix();
}

TEST(ConservativeVolumeTransferTest, transferScalarTest)
{
  auto sourceGrid = vtkUnstructuredGrid::SafeDownCast(sourceMesh->getDataSet());
  // TODO set data method within transfer
  auto pointValues = vtkSmartPointer<vtkDoubleArray>::New();
  pointValues->SetName("scalar field");
  for(int ptIdx = 0; ptIdx < sourceGrid->GetNumberOfPoints(); ++ptIdx)
  {
    double srcPt[3];
    sourceGrid->GetPoint(ptIdx, srcPt);
    // pointValues->InsertValue(ptIdx, 1.);
    pointValues->InsertValue(ptIdx, srcPt[2]);
  }
  sourceGrid->GetPointData()->AddArray(pointValues);

  std::vector<int> arrayIds = { 0 };
  volumeTransfer->transferPointData(arrayIds);

  auto end_time = std::chrono::system_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time)
                        .count();

  std::cout << "TOTAL TRANSFER TIME : " << elapsed_ms << std::endl;

  auto targetGrid = volumeTransfer->getTargetGrid();

  // TODO array indexing
  auto targetPointValues = vtkDoubleArray::SafeDownCast(targetGrid->GetPointData()->GetArray(0));
  auto errorValues = vtkSmartPointer<vtkDoubleArray>::New();
  errorValues->SetName("error");

  double maxError = 0.;

  double errorSum = 0.;
  double avgError = 0.;
  long numAccurate = 0;

  for(int ptIdx = 0; ptIdx < targetGrid->GetNumberOfPoints(); ++ptIdx)
  {
    double tgtPt[3];
    targetGrid->GetPoint(ptIdx, tgtPt);

    double expectedValue = tgtPt[2];
    double error = std::abs(targetPointValues->GetValue(ptIdx) - expectedValue);
    errorValues->InsertValue(ptIdx, error);
    errorSum += error;
    maxError = std::max(maxError, error);

    ASSERT_LT(error, 1e-14);

    if(error < 1e-14) ++numAccurate;
  }

  targetGrid->GetPointData()->AddArray(errorValues);

  avgError = errorSum/targetGrid->GetNumberOfPoints();

  std::cout << "max error : " << maxError << std::endl;
  std::cout << "avg error : " << avgError << std::endl;
  std::cout << "num accurate / total : " << numAccurate << " / " << targetGrid->GetNumberOfPoints() << std::endl;
  ASSERT_LT(maxError, 1e-13);
  ASSERT_LT(avgError, 1e-14);
}

TEST(ConservativeVolumeTransferTest, transferVectorTest)
{
  auto sourceGrid = vtkUnstructuredGrid::SafeDownCast(sourceMesh->getDataSet());

  auto pointValues = vtkSmartPointer<vtkDoubleArray>::New();
  pointValues->SetName("vector field");
  pointValues->SetNumberOfComponents(3);
  for(int ptIdx = 0; ptIdx < sourceGrid->GetNumberOfPoints(); ++ptIdx)
  {
    double srcPt[3];
    sourceGrid->GetPoint(ptIdx, srcPt);
    pointValues->InsertNextTuple3(srcPt[0], srcPt[1], srcPt[2]);
  }
  sourceGrid->GetPointData()->AddArray(pointValues);

  volumeTransfer->transfer(1);

  auto targetGrid = volumeTransfer->getTargetGrid();

  auto targetPointValues = vtkDoubleArray::SafeDownCast(targetGrid->GetPointData()->GetArray(2));
  double maxError = 0.;

  double errorSum = 0.;
  double avgError = 0.;

  for(int ptIdx = 0; ptIdx < targetGrid->GetNumberOfPoints(); ++ptIdx)
  {
    for(int componentId = 0; componentId < 3; ++componentId)
    {
      double tgtPt[3];
      targetGrid->GetPoint(ptIdx, tgtPt);
      // ASSERT_NEAR(targetPointValues->GetValue(ptIdx), 1., 1e-14);
      double expectedValue = tgtPt[componentId];
      double error = std::abs(targetPointValues->GetComponent(ptIdx, componentId) - expectedValue);
      errorSum += error;
      if(error > maxError) maxError = error;
    }
  }

  avgError = errorSum/targetGrid->GetNumberOfPoints();

  std::cout << "max error : " << maxError << std::endl;
  std::cout << "avg error : " << avgError << std::endl;
  ASSERT_LT(maxError, 1e-13);
  ASSERT_LT(avgError, 1e-14);
}

TEST(SupermeshTransferTest, sourceWriteTest)
{
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetInputData(volumeTransfer->getSourceGrid());
  writer->SetFileName("output_source.vtu");
  writer->Write();
}

TEST(ConservativeVolumeTransferTest, targetWriteTest)
{
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetInputData(volumeTransfer->getTargetGrid());
  writer->SetFileName("output_target2.vtu");
  writer->Write();
}

int main(int argc, char ** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  sourceFile = argv[1];
  targetFile = argv[2];
  outputFile = "output_target1.vtu";
  return RUN_ALL_TESTS();
}
