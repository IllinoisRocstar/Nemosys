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
  EXPECT_EQ(cmp1->getNumberOfCells(), cmp2->getNumberOfCells());

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}

TEST(rocPack, NumNodesPeriodicSpheres)
{  
  meshBase* cmp1 = meshBase::Create( "periodicGeom.vtk" );
  meshBase* cmp2 = meshBase::Create( "periodicGeom_ref.vtk" );
  EXPECT_EQ(cmp1->getNumberOfPoints(), cmp2->getNumberOfPoints());

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}

TEST(rocPack, NumCellsBoundaryPacks)
{
  auto* objrocPck = 
        new NEM::GEO::rocPack("rocOut2", "boundaryPacks");
  objrocPck->rocPack2Surf();

  if (objrocPck)
    delete objrocPck;

  meshBase* cmp1 = meshBase::Create( "boundaryPacks.vtk" );
  meshBase* cmp2 = meshBase::Create( "shapes_ref.vtk" );
  EXPECT_EQ(cmp1->getNumberOfCells(), cmp2->getNumberOfCells());

  if (cmp1)
    delete cmp1;
  if (cmp2)
    delete cmp2;
}

TEST(rocPack, NumNodesBoundaryPacks)
{  
  meshBase* cmp1 = meshBase::Create( "boundaryPacks.vtk" );
  meshBase* cmp2 = meshBase::Create( "shapes_ref.vtk" );
  EXPECT_EQ(cmp1->getNumberOfPoints(), cmp2->getNumberOfPoints());

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