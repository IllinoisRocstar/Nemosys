#include <iostream>
#include <string>
#include "blockMeshGen.H"
#include "blockMeshParams.H"
#include <boost/filesystem.hpp>
#include "meshGen.H"
#include <vtkUnstructuredGrid.h>
#include <set>
#include <cstdlib>
#include <sstream>

// OpenFOAM headers
#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"
#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"
#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"
#include "Pair.H"
#include "slidingInterface.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "vtkTopo.H"
#include "fileName.H"
#include "argList.H"
#include "IFstream.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "boundBox.H"
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

blockMeshGen::blockMeshGen() // Default constructor body
{
  // booleans
  _params = new blockMeshParams();
  defaults = true;

  // initializing the Foam environment
  initialize();
}


// Constructor with user define parameters
blockMeshGen::blockMeshGen(blockMeshParams* params):
    defaults(false), _params(params)
{
  // Initialize foam environment
  initialize();
}

blockMeshGen::~blockMeshGen() // destructor definition
{
  if (defaults)
    delete _params; // Deletes object 
}

void blockMeshGen::initialize()
{
  // creates dictionaries in constant folder
  createBlockMshDict();
  createControlDict();
  createfvSchemesDict();
  createfvSolutionDict();
}

void blockMeshGen::createControlDict()
{
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the cfMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/controlDict");
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    location  \"system\";\n\
    object    controlDict;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//\n\n\
deltaT  1;\n\n\
startTime   0;\n\n\
writeInterval   1;\n\n\
// **********************************************************************//";
  contDict << contText;
  contDict.close();
}

void blockMeshGen::createfvSchemesDict()
{
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
      boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
      std::cerr << "Problem in creating system directory for the cfMesh"
                << "\n";
      std::cerr << e.what() << std::endl;
      throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/fvSchemes");
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    object    fvSchemes;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//\n\n\
gradSchemes\n\
{\n\
    default         Gauss linear;\n\
    grad(p)         Gauss linear;\n\
}\n\
\n\
divSchemes\n\
{\n\
    default         none;\n\
    div(phi,U)      Gauss linear;\n\
}\n\
\n\
laplacianSchemes\n\
{\n\
    default         none;\n\
    laplacian(nu,U) Gauss linear corrected;\n\
    laplacian((1|A(U)),p) Gauss linear corrected;\n\
}\n\
// ************************************************************************//";
    contDict << contText;
    contDict.close();
}

void blockMeshGen::createfvSolutionDict()
{
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the cfMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/fvSolution");
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    object    fvSolution;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
\n\
// ***********************************************************************//";
  contDict << contText;
  contDict.close();
}


void blockMeshGen::createBlockMshDict()
{
  // created blockMeshDict Automatically
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the blockMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  // creating mesh dictionary file
  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/blockMeshDict");
  // header
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: blockMesh interface                   |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    location  \"system\";\n\
    object    blockMeshDict;\n\
}\n\n";

    contText = contText + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    if (auto box = std::dynamic_pointer_cast<bmBox>(_params->shape))
    {
    double finalX = 0;
    double finalY = 0;
    double finalZ = 0;

    // Automatically Generates Box Around Packs and Meshes that Hexahedrally.
    if (box->autoGenerate.has_value())
    {
      const auto &autoGenerate = box->autoGenerate.value();
      Foam::fileName inFileName((autoGenerate.packFileName));

      Foam::triSurf origSurface(inFileName);
      Foam::triSurfModifier sMod(origSurface);
      Foam::pointField& points = sMod.pointsAccess();

      const Foam::boundBox bb(points);

      Foam::vector negOffset, posOffset;

      negOffset[0] = (autoGenerate.offset[0]);
      negOffset[1] = (autoGenerate.offset[1]);
      negOffset[2] = (autoGenerate.offset[2]);
      posOffset[0] = (autoGenerate.offset[0]);
      posOffset[1] = (autoGenerate.offset[1]);
      posOffset[2] = (autoGenerate.offset[2]);

      const Foam::boundBox newBB(bb.min()-negOffset, bb.max()+posOffset);

      (box->init[0]) = newBB.min()[0];
      (box->init[1]) = newBB.min()[1];
      (box->init[2]) = newBB.min()[2];

      finalX = newBB.max()[0];
      finalY = newBB.max()[1];
      finalZ = newBB.max()[2];

      std::array<double, 3> minPointsBox{newBB.min()[0], newBB.min()[1],
                                         newBB.min()[2]};
      std::array<double, 3> maxPointsBox{newBB.max()[0], newBB.max()[1],
                                         newBB.max()[2]};

      box->coordsBox = std::make_pair(minPointsBox,maxPointsBox);

      if (_params->cellSize.has_value())
      {
        double xLength = std::sqrt(((box->init[0]) - (finalX)) *
                                   ((box->init[0]) - (finalX)));
        double yLength = std::sqrt(((box->init[1]) - (finalY)) *
                                   ((box->init[1]) - (finalY)));
        double zLength = std::sqrt(((box->init[2]) - (finalZ)) *
                                   ((box->init[2]) - (finalZ)));

        (_params->nCells[0]) = xLength/(_params->cellSize.value());
        (_params->nCells[1]) = yLength/(_params->cellSize.value());
        (_params->nCells[2]) = zLength/(_params->cellSize.value());
      }
    }
    else
    {
      finalX = (box->init[0]) + (box->len[0]);
      finalY = (box->init[1]) + (box->len[1]);
      finalZ = (box->init[2]) + (box->len[2]);

      std::array<double, 3> minPointsBox{box->init[0], box->init[1],
                                         box->init[2]};
      std::array<double, 3> maxPointsBox{finalX, finalY, finalZ};
      box->coordsBox = std::make_pair(minPointsBox, maxPointsBox);
    }

    // Box data
    contText = contText + "convertToMeters " + 
            std::to_string(_params->cnvrtToMeters) + ";\n";
    contText = contText + "\nvertices\n";
    contText = contText + "(\n\n";
    contText = contText + "\t(" + std::to_string(box->init[0])
            + " " + std::to_string(box->init[1]) + " "
            + std::to_string(box->init[2]) + ")\n";
            
    contText = contText + "\t(" + std::to_string(finalX)
            + " " + std::to_string(box->init[1]) + " "
            + std::to_string(box->init[2]) + ")\n";
            
    contText = contText + "\t(" + std::to_string(finalX)
            + " " + std::to_string(finalY) + " "
            + std::to_string(box->init[2]) + ")\n";
            
    contText = contText + "\t(" + std::to_string(box->init[0])
            + " " + std::to_string(finalY) + " "
            + std::to_string(box->init[2]) + ")\n";
            
    contText = contText + "\t(" + std::to_string(box->init[0])
            + " " + std::to_string(box->init[1]) + " "
            + std::to_string(finalZ) + ")\n";
            
    contText = contText + "\t(" + std::to_string(finalX)
            + " " + std::to_string(box->init[1]) + " "
            + std::to_string(finalZ) + ")\n";
            
    contText = contText + "\t(" + std::to_string(finalX)
            + " " + std::to_string(finalY) + " "
            + std::to_string(finalZ) + ")\n";
            
    contText = contText + "\t(" + std::to_string(box->init[0])
            + " " + std::to_string(finalY) + " "
            + std::to_string(finalZ) + ")\n";
            
    contText = contText + "\n);\n";
        
    contText = contText + "\nblocks\n(";
    contText = contText + "\n\thex (0 1 2 3 4 5 6 7) (" +
            std::to_string(_params->nCells[0]) + " " +
            std::to_string(_params->nCells[1]) + " " +
            std::to_string(_params->nCells[2]) + ") simpleGrading (" +
            std::to_string(box->smplGrading[0]) + " " +
            std::to_string(box->smplGrading[1]) + " " +
            std::to_string(box->smplGrading[2]) + ")\n\n";
        
    contText = contText + ");\n";
        
    contText = contText + "\nedges\n";
    contText = contText + "(\n\n);\n";
        
    contText = contText + "\npatches";
    contText = contText + "\n(\n";
    contText = contText + "\n\tpatch up\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(3 7 6 2)\n";
    contText = contText + "\t)\n";
        
    contText = contText + "\tpatch down\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(0 4 5 1)\n";
    contText = contText + "\t)\n";
        
    contText = contText + "\tpatch in\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(3 7 4 0)\n";
    contText = contText + "\t)\n";
        
    contText = contText + "\tpatch out\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(2 6 5 1)\n";
    contText = contText + "\t)\n";
        
    contText = contText + "\tpatch front\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(0 1 2 3)\n";
    contText = contText + "\t)\n";
        
    contText = contText + "\tpatch back\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(4 5 6 7)\n";
    contText = contText + "\t)\n";
        
    contText = contText + "\n);\n";
    
    contText = contText + "\nmergePatchPairs\n(\n);\n";
    contText = contText + 
    "\n //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//";
    }

  else if (auto sphere = std::dynamic_pointer_cast<bmSphere>(_params->shape))
  {
    // Put sphere implementation here
    contText = contText + "\ncx " + std::to_string(sphere->center[0]) + ";";
    contText = contText + "\ncy " + std::to_string(sphere->center[1]) + ";";
    contText = contText + "\ncz " + std::to_string(sphere->center[2]) + ";";
    contText = contText + "\nrad " + std::to_string(sphere->radius) + ";";

    contText = contText + "\ngeometry\n{\n";
    contText = contText + "\tsphere\n\t{\n";
    contText = contText + "\t\ttype searchableSphere;\n";
    contText = contText + "\t\tcentre ($cx $cy $cz);\n";
    contText = contText + "\t\tradius $rad;\n";
    contText = contText + "\t}\n}\n";

    contText = contText + "\nconvertToMeters "
            + std::to_string(_params->cnvrtToMeters) + ";\n";

    contText = contText + "\nvx   #calc \"$cx + 0.5773502*$rad\";\n";
    contText = contText + "mvx  #calc \"$cx + -0.5773502*$rad\";\n";

    contText = contText + "\nvy  #calc \"$cy + 0.5773502*$rad\";\n";
    contText = contText + "mvy  #calc \"$cy + -0.5773502*$rad\";\n";

    contText = contText + "\nvz  #calc \"$cz + 0.5773502*$rad\";\n";
    contText = contText + "mvz  #calc \"$cz + -0.5773502*$rad\";\n";

    contText = contText + "\nax  #calc \"$cx + 0.7071067*$rad\";\n";
    contText = contText + "max  #calc \"$cx + -0.7071067*$rad\";\n";

    contText = contText + "\nay  #calc \"$cy + 0.7071067*$rad\";\n";
    contText = contText + "may  #calc \"$cy + -0.7071067*$rad\";\n";

    contText = contText + "\naz  #calc \"$cz + 0.7071067*$rad\";\n";
    contText = contText + "maz  #calc \"$cz + -0.7071067*$rad\";\n";


    contText = contText + "\nvertices\n";
    contText = contText + "(\n";
    contText = contText + "\t($mvx $mvy $mvz)\n";
    contText = contText + "\t( $vx $mvy $mvz)\n";
    contText = contText + "\t( $vx  $vy $mvz)\n";
    contText = contText + "\t($mvx  $vy $mvz)\n";
    contText = contText + "\t($mvx $mvy  $vz)\n";
    contText = contText + "\t( $vx $mvy  $vz)\n";
    contText = contText + "\t( $vx  $vy  $vz)\n";
    contText = contText + "\t($mvx  $vy  $vz)\n";
    contText = contText + ");\n";


    contText = contText + "\nblocks\n";
    contText = contText + "(\n";
    contText = contText + "\thex (0 1 2 3 4 5 6 7) (" +
            std::to_string(_params->nCells[0]) + " " +
            std::to_string(_params->nCells[1]) + " " +
            std::to_string(_params->nCells[2]) + ") simpleGrading (" +
            std::to_string(sphere->sphrGrading[0]) + " " +
            std::to_string(sphere->sphrGrading[1]) + " " +
            std::to_string(sphere->sphrGrading[2]) + ")\n\n";

    contText = contText + ");\n";


    contText = contText + "\nedges\n";
    contText = contText + "(\n";
    contText = contText + "\tarc 0 1 ($cx $may $maz)\n";
    contText = contText + "\tarc 2 3 ($cx $ay  $maz)\n";
    contText = contText + "\tarc 6 7 ($cx $ay  $az)\n";
    contText = contText + "\tarc 4 5 ($cx $may $az)\n";

    contText = contText + "\tarc 0 3 ($max $cy $maz)\n";
    contText = contText + "\tarc 1 2 ($ax  $cy $maz)\n";
    contText = contText + "\tarc 5 6 ($ax  $cy $az)\n";
    contText = contText + "\tarc 4 7 ($max $cy $az)\n";

    contText = contText + "\tarc 0 4 ($max $may $cz)\n";
    contText = contText + "\tarc 1 5 ($ax  $may $cz)\n";
    contText = contText + "\tarc 2 6 ($ax  $ay  $cz)\n";
    contText = contText + "\tarc 3 7 ($max $ay  $cz)\n";

    contText = contText + ");\n";

    contText = contText + "\nfaces\n(\n";
    contText = contText + "\tproject (0 4 7 3) sphere\n";
    contText = contText + "\tproject (2 6 5 1) sphere\n";
    contText = contText + "\tproject (1 5 4 0) sphere\n";
    contText = contText + "\tproject (3 7 6 2) sphere\n";
    contText = contText + "\tproject (0 3 2 1) sphere\n";
    contText = contText + "\tproject (4 5 6 7) sphere\n";
    contText = contText + ");\n";

    contText = contText + "\nboundary\n(\n";
    contText = contText + "\twalls\n\t{\n";
    contText = contText + "\t\ttype wall;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(0 4 7 3)\n";
    contText = contText + "\t\t\t(2 6 5 1)\n";
    contText = contText + "\t\t\t(1 5 4 0)\n";
    contText = contText + "\t\t\t(3 7 6 2)\n";
    contText = contText + "\t\t\t(0 3 2 1)\n";
    contText = contText + "\t\t\t(4 5 6 7)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";
    contText = contText + ");\n";

    contText = contText + 
    "\n //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//";
  }
  else if (auto cylTaperedCone
             = std::dynamic_pointer_cast<bmCylTaperedCone>(_params->shape))
  {
    contText = contText + "\n\ncx " +
               std::to_string(cylTaperedCone->centerCyl[0]) + ";\n";
    contText =
        contText + "cy " + std::to_string(cylTaperedCone->centerCyl[1]) + ";\n";
    contText =
        contText + "cz " + std::to_string(cylTaperedCone->centerCyl[2]) + ";\n";
    contText =
        contText + "rad1 " + std::to_string(cylTaperedCone->radius1) + ";\n";
    contText = contText + "rad2 " +
               std::to_string(
                   cylTaperedCone->radius2.value_or(cylTaperedCone->radius1)) +
               ";\n";
    contText = contText + "h " + std::to_string(cylTaperedCone->height) + ";\n";
    contText = contText + "cellX " + std::to_string(_params->nCells[0]) + ";\n";
    contText = contText + "cellY " + std::to_string(_params->nCells[1]) + ";\n";
    contText = contText + "cellZ " + std::to_string(_params->nCells[2]) + ";\n";
    contText = contText + "grdX " +
               std::to_string(cylTaperedCone->cylGrading[0]) + ";\n";
    contText = contText + "grdY " +
               std::to_string(cylTaperedCone->cylGrading[1]) + ";\n";
    contText = contText + "grdZ " +
               std::to_string(cylTaperedCone->cylGrading[2]) + ";\n\n";

    contText = contText + "\nX0 #calc \"$cx - 0.707106*$rad1\";";
    contText = contText + "\nY0 #calc \"$cy - 0.707106*$rad1\";";
    contText = contText + "\nX1 #calc \"$cx + 0.707106*$rad1\";";
    contText = contText + "\nY1 #calc \"$cy - 0.707106*$rad1\";";
    contText = contText + "\nX2 #calc \"$cx - 0.707106*0.3*$rad1\";";
    contText = contText + "\nY2 #calc \"$cy - 0.707106*0.3*$rad1\";";
    contText = contText + "\nX3 #calc \"$cx + 0.707106*0.3*$rad1\";";
    contText = contText + "\nY3 #calc \"$cy - 0.707106*0.3*$rad1\";";
    contText = contText + "\nX4 #calc \"$cx + 0.707106*0.3*$rad1\";";
    contText = contText + "\nY4 #calc \"$cy + 0.707106*0.3*$rad1\";";
    contText = contText + "\nX5 #calc \"$cx + 0.707106*$rad1\";";
    contText = contText + "\nY5 #calc \"$cy + 0.707106*$rad1\";";
    contText = contText + "\nX6 #calc \"$cx - 0.707106*0.3*$rad1\";";
    contText = contText + "\nY6 #calc \"$cy + 0.707106*0.3*$rad1\";";
    contText = contText + "\nX7 #calc \"$cx - 0.707106*$rad1\";";
    contText = contText + "\nY7 #calc \"$cy + 0.707106*$rad1\";";
    contText = contText + "\nX8 #calc \"$cx - 0.707106*$rad2\";";
    contText = contText + "\nY8 #calc \"$cy - 0.707106*$rad2\";";
    contText = contText + "\nX9 #calc \"$cx + 0.707106*$rad2\";";
    contText = contText + "\nY9 #calc \"$cy - 0.707106*$rad2\";";
    contText = contText + "\nX10 #calc \"$cx - 0.707106*0.3*$rad2\";";
    contText = contText + "\nY10 #calc \"$cy - 0.707106*0.3*$rad2\";";
    contText = contText + "\nX11 #calc \"$cx + 0.707106*0.3*$rad2\";";
    contText = contText + "\nY11 #calc \"$cy - 0.707106*0.3*$rad2\";";
    contText = contText + "\nX12 #calc \"$cx + 0.707106*0.3*$rad2\";";
    contText = contText + "\nY12 #calc \"$cy + 0.707106*0.3*$rad2\";";
    contText = contText + "\nX13 #calc \"$cx + 0.707106*$rad2\";";
    contText = contText + "\nY13 #calc \"$cy + 0.707106*$rad2\";";
    contText = contText + "\nX14 #calc \"$cx - 0.707106*0.3*$rad2\";";
    contText = contText + "\nY14 #calc \"$cy + 0.707106*0.3*$rad2\";";
    contText = contText + "\nX15 #calc \"$cx - 0.707106*$rad2\";";
    contText = contText + "\nY15 #calc \"$cy + 0.707106*$rad2\";";
    contText = contText + "\narcX1 #calc \"$cx*1\";";
    contText = contText + "\narcY1 #calc \"$cy - $rad1\";";
    contText = contText + "\narcX2 #calc \"$cx + $rad1\";";
    contText = contText + "\narcY2 #calc \"$cy*1\";";
    contText = contText + "\narcX3 #calc \"$cx*1\";";
    contText = contText + "\narcY3 #calc \"$cy + $rad1\";";
    contText = contText + "\narcX4 #calc \"$cx - $rad1\";";
    contText = contText + "\narcY4 #calc \"$cy*1\";";
    contText = contText + "\narcX5 #calc \"$cx*1\";";
    contText = contText + "\narcY5 #calc \"$cy - $rad2\";";
    contText = contText + "\narcX6 #calc \"$cx + $rad2\";";
    contText = contText + "\narcY6 #calc \"$cy*1\";";
    contText = contText + "\narcX7 #calc \"$cx*1\";";
    contText = contText + "\narcY7 #calc \"$cy + $rad2\";";
    contText = contText + "\narcX8 #calc \"$cx - $rad2\";";
    contText = contText + "\narcY8 #calc \"$cy*1\";";
    contText = contText + "\nfz #calc \"$cz + $h\";";
    contText = contText + "\nconvertToMeters "
            + std::to_string(_params->cnvrtToMeters) + ";\n";
    
    contText = contText + "\nvertices\n(\n\n";
    contText = contText + "\t($X0 $Y0 $cz)\n";
    contText = contText + "\t($X1 $Y1 $cz)\n";
    contText = contText + "\t($X2 $Y2 $cz)\n";
    contText = contText + "\t($X3 $Y3 $cz)\n";
    contText = contText + "\t($X4 $X4 $cz)\n";
    contText = contText + "\t($X5 $Y5 $cz)\n";
    contText = contText + "\t($X6 $Y6 $cz)\n";
    contText = contText + "\t($X7 $Y7 $cz)\n";
    contText = contText + "\t($X8 $Y8 $fz)\n";
    contText = contText + "\t($X9 $Y9 $fz)\n";
    contText = contText + "\t($X10 $Y10 $fz)\n";
    contText = contText + "\t($X11 $Y11 $fz)\n";
    contText = contText + "\t($X12 $X12 $fz)\n";
    contText = contText + "\t($X13 $Y13 $fz)\n";
    contText = contText + "\t($X14 $Y14 $fz)\n";
    contText = contText + "\t($X15 $Y15 $fz)\n";
    contText = contText + ");\n\n";
    
    contText = contText
            + "\nblocks\n(\n";
    contText = contText
            + "\thex (2 3 4 6 10 11 12 14) ($cellX $cellY $cellZ) "
            +"simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText
            + "\thex (0 1 3 2 8 9 11 10) ($cellX $cellY $cellZ) "
            +"simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText
            + "\thex (3 1 5 4 11 9 13 12) ($cellX $cellY $cellZ) "
            +"simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText
            + "\thex (4 5 7 6 12 13 15 14) ($cellX $cellY $cellZ) "
            +"simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText
            + "\thex (6 7 0 2 14 15 8 10) ($cellX $cellY $cellZ) "
            +"simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText
            + ");\n";
    
    contText = contText + "\nedges\n(\n";
    contText = contText + "\tarc 0 1 ($arcX1 $arcY1 $cz)\n";
    contText = contText + "\tarc 1 5 ($arcX2 $arcY2 $cz)\n";
    contText = contText + "\tarc 5 7 ($arcX3 $arcY3 $cz)\n";
    contText = contText + "\tarc 7 0 ($arcX4 $arcY4 $cz)\n";
    contText = contText + "\tarc 8 9 ($arcX5 $arcY5 $fz)\n";
    contText = contText + "\tarc 9 13 ($arcX6 $arcY6 $fz)\n";
    contText = contText + "\tarc 13 15 ($arcX7 $arcY7 $fz)\n";
    contText = contText + "\tarc 15 8 ($arcX8 $arcY8 $fz)\n";
    contText = contText + ");\n\n";
    
    contText = contText + "boundary\n(\n";
    contText = contText + "\tinlet\n\t{\n";
    contText = contText + "\t\ttype patch;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(0 1 2 3)\n";
    contText = contText + "\t\t\t(1 3 4 5)\n";
    contText = contText + "\t\t\t(4 5 7 6)\n";
    contText = contText + "\t\t\t(0 2 6 7)\n";
    contText = contText + "\t\t\t(2 3 4 6)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";
    
    contText = contText + "\toutlet\n\t{\n";
    contText = contText + "\t\ttype patch;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(8 9 10 11)\n";
    contText = contText + "\t\t\t(11 9 12 13)\n";
    contText = contText + "\t\t\t(12 13 14 15)\n";
    contText = contText + "\t\t\t(14 15 10 8)\n";
    contText = contText + "\t\t\t(10 11 14 12)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";
    
    contText = contText + "\twalls\n\t{\n";
    contText = contText + "\t\ttype patch;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(1 5 13 9)\n";
    contText = contText + "\t\t\t(5 7 15 13)\n";
    contText = contText + "\t\t\t(7 0 8 15)\n";
    contText = contText + "\t\t\t(0 1 9 8)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";
    contText = contText + ");\n";
    
    contText = contText + 
        "\n //* * * * * * * * * * * * * * * * * * * * * * * * * * * * *//";
  }
  else
  {
    std::cerr << "User cannot select multiple geometries in single run!\n" 
              << std::endl;
    exit(1);
  }
    
  contDict << contText;   
  contDict.close();
  
}

// Implementation of blockMesh code
int blockMeshGen::createMeshFromSTL(const char* fname)
{

  using namespace Foam;

  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info<< "Create time\n" << Foam::endl;
  Foam::argList::noParallel();

  Time runTime
  (
    Time::controlDictName,
    "",
    ""
  );


  word regionName;
  word regionPath;
  regionName = polyMesh::defaultRegion;
  regionPath = regionName;
  
  // Locating blockMeshDict in system directory
  const word dictName("blockMeshDict");
  fileName dictPath;
  dictPath = runTime.system()/dictName;
  
  // Looks for and cleans polyMesh (Can use boost directory if needed)
  fileName polyMeshPath
  (
    runTime.path()/runTime.constant()/polyMesh::meshSubDir
  );

  if (exists(polyMeshPath))
  {
    if (exists(polyMeshPath/dictName))
    {
      Info<< "Not deleting polyMesh directory " << nl
          << "    " << polyMeshPath << nl
          << "    because it contains " << dictName << endl;
    }
    else
    {
      Info<< "Deleting polyMesh directory" << nl
          << "    " << polyMeshPath << endl;
      rmDir(polyMeshPath);
    }
  } 

  // Creates blockMesh from defined dictionary
  IOobject meshDictIO
  (
    dictPath,
    runTime,
    IOobject::MUST_READ,
    IOobject::NO_WRITE,
    false
  );

#ifdef HAVE_OF5

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true)) //OF-5.0
    //if (!meshDictIO.headerOk())  // OF-4.0
    {
        FatalErrorInFunction
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

      Info<< "Creating block mesh from\n    "
      << meshDictIO.objectPath() << endl;

  IOdictionary meshDict(meshDictIO);
  blockMesh blocks(meshDict, regionName);
  
  Info<< nl << "Creating polyMesh from blockMesh" << endl;

  word defaultFacesName = "defaultFaces";
  word defaultFacesType = emptyPolyPatch::typeName;
  polyMesh mesh
  (
    IOobject
    (
      regionName,
      runTime.constant(),
      runTime
    ),
    xferCopy(blocks.points()),
    blocks.cells(),
    blocks.patches(),
    blocks.patchNames(),
    blocks.patchDicts(),
    defaultFacesName,
    defaultFacesType
  );

#endif

#ifdef HAVE_OF4

    //if (!meshDictIO.typeHeaderOk<IOdictionary>(true)) //OF-5.0
    if (!meshDictIO.headerOk())  // OF-4.0
    {
        FatalErrorInFunction
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

      Info<< "Creating block mesh from\n    "
      << meshDictIO.objectPath() << endl;

  IOdictionary meshDict(meshDictIO);
  blockMesh blocks(meshDict, regionName);
  
  Info<< nl << "Creating polyMesh from blockMesh" << endl;

  word defaultFacesName = "defaultFaces";
  word defaultFacesType = emptyPolyPatch::typeName;
  polyMesh mesh
  (
    IOobject
    (
      regionName,
      runTime.constant(),
      runTime
    ),
    xferCopy(blocks.points()),
    blocks.cells(),
    blocks.patches(),
    blocks.patchNames(),
    blocks.patchDicts(),
    defaultFacesName,
    defaultFacesType
  );
    
#endif

#ifdef HAVE_OF6

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true)) //OF-5.0
    //if (!meshDictIO.headerOk())  // OF-4.0
    {
        FatalErrorInFunction
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

      Info<< "Creating block mesh from\n    "
      << meshDictIO.objectPath() << endl;

  IOdictionary meshDict(meshDictIO);
  blockMesh blocks(meshDict, regionName);
  
  Info<< nl << "Creating polyMesh from blockMesh" << endl;

  word defaultFacesName = "defaultFaces";
  word defaultFacesType = emptyPolyPatch::typeName;
  polyMesh mesh
  (
    IOobject
    (
      regionName,
      runTime.constant(),
      runTime
    ),
    xferCopy(blocks.points()),
    blocks.cells(),
    blocks.patches(),
    blocks.patchNames(),
    blocks.patchDicts(),
    defaultFacesName,
    defaultFacesType
  );

#endif

#ifdef HAVE_OF7

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true)) //OF-5.0
    //if (!meshDictIO.headerOk())  // OF-4.0
    {
        FatalErrorInFunction
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from\n    "
    << meshDictIO.objectPath() << endl;

    IOdictionary meshDict(meshDictIO);
    blockMesh blocks(meshDict, regionName);
  
    Info<< nl << "Creating polyMesh from blockMesh" << endl;

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    polyMesh mesh
    (
        IOobject
        (
        regionName,
        runTime.constant(),
        runTime
        ),
        Foam::clone(blocks.points()),
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        defaultFacesName,
        defaultFacesType
    );

#endif

    // Disabling mergePairs feature for now
    /*// Read in a list of dictionaries for the merge patch pairs
    if (meshDict.found("mergePatchPairs"))
    {
        List<Pair<word>> mergePatchPairs
        (
            meshDict.lookup("mergePatchPairs")
        );

        #include "mergePatchPairs.H"
    }
    else
    {
        Info<< nl << "There are no merge patch pairs edges" << endl;
    }*/
  
  label nZones = blocks.numZonedBlocks();

  if (nZones > 0)
  {
    Info<< nl << "Adding cell zones" << endl;

    // Map from zoneName to cellZone index
    HashTable<label> zoneMap(nZones);

    // Cells per zone.
    List<DynamicList<label>> zoneCells(nZones);

    // Running cell counter
    label celli = 0;

    // Largest zone so far
    label freeZoneI = 0;

    forAll(blocks, blockI)
    {
      const block& b = blocks[blockI];

      #ifdef HAVE_OF5

        const List<FixedList<label, 8>> blockCells = b.cells();  //OF-5.0

      #endif
      #ifdef HAVE_OF4

        const labelListList& blockCells = b.cells();         //OF-4.0

      #endif
      #ifdef HAVE_OF6

        const List<FixedList<label, 8>> blockCells = b.cells();  //OF-5.0

      #endif
      #ifdef HAVE_OF7

        const List<FixedList<label, 8>> blockCells = b.cells();  //OF-5.0

      #endif

       const word& zoneName = b.zoneName();

       if (zoneName.size())
       {
          HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

          label zoneI;

          if (iter == zoneMap.end())
          {
            zoneI = freeZoneI++;

            Info<< "    " << zoneI << '\t' << zoneName << endl;

            zoneMap.insert(zoneName, zoneI);
          }
          else
          {
            zoneI = iter();
          }

          forAll(blockCells, i)
          {
            zoneCells[zoneI].append(celli++);
          }

        }
        else
        {
          celli += b.cells().size();
        }

    }


    List<cellZone*> cz(zoneMap.size());

    Info<< nl << "Writing cell zones as cellSets" << endl;

    forAllConstIter(HashTable<label>, zoneMap, iter)
    {
      label zoneI = iter();

      cz[zoneI] = new cellZone
      (
        iter.key(),
        zoneCells[zoneI].shrink(),
        zoneI,
        mesh.cellZones()
      );

      // Write as cellSet for ease of processing
      cellSet cset(mesh, iter.key(), zoneCells[zoneI].shrink());
      cset.write();
    }

    mesh.pointZones().setSize(0);
    mesh.faceZones().setSize(0);
    mesh.cellZones().setSize(0);
    mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
  }

  // Set the precision of the points data to 10
  IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

  Info<< nl << "Writing polyMesh" << endl;
  mesh.removeFiles();
  if (!mesh.write())
  {
    FatalErrorInFunction
        << "Failed writing polyMesh."
        << exit(FatalError);
  }


  //
  // writes some information
  //
  {
    const polyPatchList& patches = mesh.boundaryMesh();

    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
        << "  " << "nPoints: " << mesh.nPoints() << nl
        << "  " << "nCells: " << mesh.nCells() << nl
        << "  " << "nFaces: " << mesh.nFaces() << nl
        << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

    Info<< "----------------" << nl
        << "Patches" << nl
        << "----------------" << nl;

    forAll(patches, patchi)
    {
      const polyPatch& p = patches[patchi];

      Info<< "  " << "patch " << patchi
          << " (start: " << p.start()
          << " size: " << p.size()
          << ") name: " << p.name()
          << nl;
    }

  }

  Info<< "\nEnd\n" << endl;

  readFoamMesh();
}

// Reads mesh from polyMesh
void blockMeshGen::readFoamMesh()
{ 
  Foam::Info<< "Create time\n" << Foam::endl;
  Foam::argList::noParallel();

  Time runTime
  (
    Time::controlDictName,
    "",
    ""
  );
  // reading mesh database and converting
  Foam::word regionName;
  regionName = Foam::fvMesh::defaultRegion;
  Foam::Info
      << "Create mesh for time = "
      << runTime.timeName() << Foam::nl << Foam::endl;

  _fmesh = new Foam::fvMesh 
  (
    Foam::IOobject
    (
      regionName,
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
  );  

  genMshDB();
}

// Generates VTK database.
void blockMeshGen::genMshDB()
{
  // Obtaining the mesh data and generating requested output file
  std::cout << "Number of points " << _fmesh->nPoints() << std::endl;
  std::cout << "Number of cells "<< _fmesh->nCells() << std::endl;

  // declare vtk dataset
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp 
      = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // decomposition
  Foam::vtkTopo::decomposePoly = false;

  // creating equivalent vtk topology from fvMesh
  // by default polyhedral cells will be decomposed to 
  // tets and pyramids. Additional points will be added
  // to underlying fvMesh.
  std::cout << "Performing topological decomposition.\n";
  Foam::vtkTopo topo(*_fmesh);

    // point coordinates
  Foam::pointField pf = _fmesh->points();
  vtkSmartPointer<vtkPoints> points 
      = vtkSmartPointer<vtkPoints>::New(); 
  for (int ipt=0; ipt<_fmesh->nPoints(); ipt++)
      points->InsertNextPoint(
              pf[ipt].x(), 
              pf[ipt].y(), 
              pf[ipt].z() 
              );
  dataSet_tmp->SetPoints(points);

  // cell types
  std::vector<int> pntIds;
  int nCelPnts = 0;
  for (int icl=0; icl<topo.vertLabels().size(); icl++)
  {
    if ( topo.cellTypes()[icl] != VTK_POLYHEDRON )
    {
      nCelPnts = topo.vertLabels()[icl].size();
      pntIds.resize(nCelPnts, -1);
      for (int ip=0; ip< nCelPnts; ip++)
        pntIds[ip] = topo.vertLabels()[icl][ip];
      createVtkCell(dataSet_tmp, topo.cellTypes()[icl], pntIds);
    }
    else
    {
      // polyhedral cells treated differently in vtk
      // faces should be defined for them
      int nFace = topo.vertLabels()[icl][0]; 
      vtkSmartPointer<vtkCellArray> faces =
          vtkSmartPointer<vtkCellArray>::New();
      std::vector<vtkIdType> faceIds;
      std::set<vtkIdType> pntIds;
      int dataId = 1;
      for (int iFace=0; iFace<nFace; iFace++)
      {
        faceIds.clear();
        int nFaceId = topo.vertLabels()[icl][dataId++];
        int pntId;
        for (int ifid=0; ifid<nFaceId; ifid++)
        {
          pntId = topo.vertLabels()[icl][dataId++];
          pntIds.insert(pntId);
          faceIds.push_back(pntId);
        }
        faces->InsertNextCell(nFaceId, &faceIds[0]);
      }
      std::vector<vtkIdType> pntIdsVec(pntIds.begin(), pntIds.end());
      dataSet_tmp->InsertNextCell(VTK_POLYHEDRON, pntIdsVec.size(),
        &pntIdsVec[0], nFace, faces->GetPointer());
    }
  }
  dataSet = dataSet_tmp;
}


void blockMeshGen::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                              const int cellType, 
                              std::vector<int>& vrtIds)
{
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(vrtIds.size());
    for (auto pit = vrtIds.begin(); pit!= vrtIds.end(); pit++)
        vtkCellIds->SetId(pit-vrtIds.begin(), *pit);
    dataSet->InsertNextCell(cellType,vtkCellIds);
}
