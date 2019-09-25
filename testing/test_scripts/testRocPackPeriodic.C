#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <gtest.h>
#include "rocPack.H"
#include "meshBase.H"


TEST(rocPack, NumCellsPeriodicSpheres)
{
  auto* objrocPck = new NEM::GEO::rocPack("rocOut", "periodicGeom");
  objrocPck->rocPack2Surf();

  if (objrocPck)
    delete objrocPck;

  meshBase* cmp1 = meshBase::Create( "periodicGeom.vtk" );
  meshBase* cmp2 = meshBase::Create( "periodicGeom_ref.vtk" );
  EXPECT_TRUE((cmp1->getNumberOfCells() >= cmp2->getNumberOfCells()*0.90) &&
              (cmp1->getNumberOfCells() <= cmp2->getNumberOfCells()*1.1));

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}

TEST(rocPack, NumNodesPeriodicSpheres)
{  
  meshBase* cmp1 = meshBase::Create( "periodicGeom.vtk" );
  meshBase* cmp2 = meshBase::Create( "periodicGeom_ref.vtk" );
  EXPECT_TRUE((cmp1->getNumberOfPoints() >= cmp2->getNumberOfPoints()*0.90) &&
              (cmp1->getNumberOfPoints() <= cmp2->getNumberOfPoints()*1.1));

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}

TEST(rocPack, NumCellsBoundaryPacks)
{
  auto* objrocPck = 
        new NEM::GEO::rocPack("rocOut", "boundaryPacks");
  objrocPck->rocPack2Surf();

  if (objrocPck)
    delete objrocPck;

  meshBase* cmp1 = meshBase::Create( "boundaryPacks.vtk" );
  meshBase* cmp2 = meshBase::Create( "shapes_ref.vtk" );
  EXPECT_TRUE((cmp1->getNumberOfCells() >= cmp2->getNumberOfCells()*0.90) &&
              (cmp1->getNumberOfCells() <= cmp2->getNumberOfCells()*1.1));

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}

TEST(rocPack, NumNodesBoundaryPacks)
{  
  meshBase* cmp1 = meshBase::Create( "boundaryPacks.vtk" );
  meshBase* cmp2 = meshBase::Create( "shapes_ref.vtk" );
  EXPECT_TRUE((cmp1->getNumberOfPoints() >= cmp2->getNumberOfPoints()*0.90) &&
              (cmp1->getNumberOfPoints() <= cmp2->getNumberOfPoints()*1.1));

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}


// test constructor
int main(int argc, char** argv) 
{
  // IO
  ::testing::InitGoogleTest(&argc, argv);

  // running tests
  int res = RUN_ALL_TESTS();

  return res;
}