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
#include <gtest/gtest.h>
#include <chrono>

#include <PatchRecovery/patchRecovery.H>
#include <Mesh/meshBase.H>

const char* nodeMesh;
const char* recoveredMesh;
const char* hexmesh;

TEST(PatchRecoveryConstructor, ConstructWithArray)
{
  // load reference node mesh
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(nodeMesh);
  const std::vector<int> arrayIDs = {0,1,2,3,4,5,7};
  int order = 1;
  std::unique_ptr<PatchRecovery> recoverObj
    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh->getDataSet(),order,arrayIDs));
  recoverObj->computeNodalError();
  mesh->write("errorTest.vtu");
}

//TEST(PatchRecoveryTensorProdBasis, RecoverNodalSol)
//{
//  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(hexmesh);
//  const std::vector<int> arrayIDs = {0,1};
//  int order = 1;
//  std::unique_ptr<PatchRecovery> recoverObj
//    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh.get(),order,arrayIDs));
//  Timer T;
//  T.start();
//  recoverObj->recoverNodalSolution(0);
//  T.stop();
//  std::cout << "time for regpoly (ms): " << T.elapsed() << std::endl;
//  mesh->write("polyApproxRecovery.vtu");
//  T.start();
//  recoverObj->recoverNodalSolution(1);
//  T.stop();
//  std::cout << "time for orthopoly(ms): " << T.elapsed() << std::endl;
//  mesh->write("orthoPolyApproxRecovery.vtu");
//}
 
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 4);
  nodeMesh = argv[1];
  recoveredMesh = argv[2];
  hexmesh = argv[3];
  return RUN_ALL_TESTS();
}
