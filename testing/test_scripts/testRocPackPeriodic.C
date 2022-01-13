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
#ifdef _WIN32
#  define _USE_MATH_DEFINES
#endif
#include <gtest/gtest.h>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <Mesh/meshBase.H>
#include <Geometry/rocPack.H>

TEST(rocPack, NumCellsPeriodicSpheres) {
  auto *objrocPck = new NEM::GEO::rocPack("rocOut", "periodicGeom.vtu");
  objrocPck->rocPack2Surf();

  if (objrocPck) delete objrocPck;

  meshBase *cmp1 = meshBase::Create("periodicGeom_oldMSH.vtu");
  meshBase *cmp2 = meshBase::Create("periodicGeom_ref.vtu");
  EXPECT_TRUE((cmp1->getNumberOfCells() >= cmp2->getNumberOfCells() * 0.90) &&
              (cmp1->getNumberOfCells() <= cmp2->getNumberOfCells() * 1.1));

  if (cmp1) delete cmp1;
  if (cmp2) delete cmp2;
}

TEST(rocPack, NumNodesPeriodicSpheres) {
  meshBase *cmp1 = meshBase::Create("periodicGeom_oldMSH.vtu");
  meshBase *cmp2 = meshBase::Create("periodicGeom_ref.vtu");
  EXPECT_TRUE((cmp1->getNumberOfPoints() >= cmp2->getNumberOfPoints() * 0.90) &&
              (cmp1->getNumberOfPoints() <= cmp2->getNumberOfPoints() * 1.1));

  if (cmp1) delete cmp1;
  if (cmp2) delete cmp2;
}

TEST(rocPack, NumCellsBoundaryPacks) {
  auto *objrocPck = new NEM::GEO::rocPack("rocOut2", "shapes.vtu");
  objrocPck->rocPack2Surf();

  if (objrocPck) delete objrocPck;

  meshBase *cmp1 = meshBase::Create("shapes_oldMSH.vtu");
  meshBase *cmp2 = meshBase::Create("shapes_ref.vtu");
  EXPECT_TRUE((cmp1->getNumberOfCells() >= cmp2->getNumberOfCells() * 0.90) &&
              (cmp1->getNumberOfCells() <= cmp2->getNumberOfCells() * 1.1));

  if (cmp1) delete cmp1;
  if (cmp2) delete cmp2;
}

TEST(rocPack, NumNodesBoundaryPacks) {
  meshBase *cmp1 = meshBase::Create("shapes_oldMSH.vtu");
  meshBase *cmp2 = meshBase::Create("shapes_ref.vtu");
  EXPECT_TRUE((cmp1->getNumberOfPoints() >= cmp2->getNumberOfPoints() * 0.90) &&
              (cmp1->getNumberOfPoints() <= cmp2->getNumberOfPoints() * 1.1));

  if (cmp1) delete cmp1;
  if (cmp2) delete cmp2;
}

TEST(rocPack, PeriodicMesh) {
  auto *objrocPck = new NEM::GEO::rocPack("rocOut_Mesh", "meshPeriodic3D.vtu");
  objrocPck->setPeriodicGeometry();
  objrocPck->sanityCheckOn();
  objrocPck->enableCohesiveElements();
  objrocPck->shrinkVolumes(0.8);
  objrocPck->setMeshSize(0.08);
  objrocPck->enableTwoPhysGrps();
  objrocPck->setPeriodicMesh();
  objrocPck->rocPack2Periodic3D();
  objrocPck->setNodeLocations(15, 30, 55);
  objrocPck->setRandomSurface(25);

  EXPECT_EQ(true, objrocPck->getTestResult());

  if (objrocPck) delete objrocPck;
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);

  // running tests
  int res = RUN_ALL_TESTS();

  return res;
}
