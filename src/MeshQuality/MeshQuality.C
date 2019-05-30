#include <MeshQuality.H>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCell.h>
#include <vtkCellType.h>
#include <AuxiliaryFunctions.H>

#ifdef HAVE_CFMSH 
// openfoam headers
#include "fvCFD.H"
#include "fvMesh.H"
#include "fileName.H"

// cfmesh headers
#include "polyMeshGenModifier.H"
#include "meshOptimizer.H"
#endif


MeshQuality::MeshQuality(const meshBase* _mesh)
  : mesh(_mesh)
{
  qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
  qualityFilter->SetInputData(mesh->getDataSet());
  qualityFilter->SetTriangleQualityMeasureToShape();
  qualityFilter->SetTetQualityMeasureToShape();
  qualityFilter->SetQuadQualityMeasureToShape();
  qualityFilter->SetHexQualityMeasureToShape();
  qualityFilter->Update();
}

MeshQuality::~MeshQuality()
{
//  mesh->unsetCellDataArray("Quality");
//  mesh->unsetFieldDataArray("Mesh Triangle Quality");
//  mesh->unsetFieldDataArray("Mesh Quadrilateral Quality");
//  mesh->unsetFieldDataArray("Mesh Tetrahedron Quality");
//  mesh->unsetFieldDataArray("Mesh Hexahedron Quality");
}

void MeshQuality::checkMesh(std::ostream& outputStream)
{

  
  outputStream << "------------- Shape Quality Statistics -------------\n" << std::endl;
  outputStream << "Cell Type" << std::setw(16) << "Num Cells" << std::setw(16)
               << "Minimum" << std::setw(16) << "Maximum" << std::setw(16)
               << "Average" << std::setw(16) << "Variance" << std::endl;
  
  for (int i = 0; i < 4; ++i)
  {
    if (i == 0) 
      outputStream << "Tri"  << std::setw(16);
    else if (i == 1) 
      outputStream << "Quad" << std::setw(15);
    else if (i == 2) 
      outputStream << "Tet" << std::setw(16);
    else 
      outputStream << "Hex" << std::setw(16);
    
    vtkSmartPointer<vtkDoubleArray> qualityField = getStats(i);
    double* val = qualityField->GetTuple(0);

    outputStream << std::right << val[4] << std::setw(16)
                 << std::right << val[0] << std::setw(16)
                 << std::right << val[2] << std::setw(16)
                 << std::right << val[1] << std::setw(16)
                 << std::right << val[3] << std::endl;

  }

  outputStream << "\n------------- Detailed Statistics ------------------\n" << std::endl;



  vtkSmartPointer<vtkDoubleArray> qualityArray =
    vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));

  mesh->getDataSet()->GetCellData()->AddArray(qualityArray);
  std::string qfn = trim_fname(mesh->getFileName(),"") + "-qal.vtu";
  mesh->write(qfn);


  // writing to file

  
  outputStream << "Type" << std::setw(10) << "Quality\n" << std::endl;
  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
  {
    int cellType = mesh->getDataSet()->GetCell(i)->GetCellType();
    double val = qualityArray->GetValue(i);
    switch(cellType)
    {
      case VTK_TRIANGLE:
      {
        outputStream << "TRI\t" << std::right << setw(10) << val << std::endl; 
        break;
      }
      case VTK_TETRA:
      {
        outputStream << "TET\t" << std::right << setw(10) << val << std::endl; 
        break;
      }
      case VTK_QUAD:
      {
        outputStream << "QUAD\t" << std::right << setw(10) << val << std::endl; 
        break;
      }
      case VTK_HEXAHEDRON:
      {
        outputStream << "HEX\t" << std::right << setw(10) << val << std::endl; 
        break;
      }
      default:
      {
        std::cout << "Invalid cell type. Only tri,quad,tet and hex supported." << std::endl;
        exit(1);
      }
      
    }
  }
}

void MeshQuality::checkMesh()
{
  checkMesh(std::cout);
}

void MeshQuality::checkMesh(const std::string &fname)
{
  std::ofstream outputStream(fname);
  if (!outputStream.good())
  {
    std::cout << "Error opening file " << fname << std::endl;
    exit(1);
  }
  checkMesh(outputStream);
}

vtkSmartPointer<vtkDoubleArray> MeshQuality::getStats(int n)
{
  vtkSmartPointer<vtkDoubleArray> qualityField = vtkSmartPointer<vtkDoubleArray>::New();
  switch(n)
  {
    case(0):
    {
      qualityField = 
        vtkDoubleArray::SafeDownCast(
          qualityFilter->GetOutput()->GetFieldData()->GetArray("Mesh Triangle Quality"));
      break;
    }
    case(1):
    {
      qualityField = 
        vtkDoubleArray::SafeDownCast(
          qualityFilter->GetOutput()->GetFieldData()->GetArray("Mesh Quadrilateral Quality"));
      break;
    }
    case(2):
    {
      qualityField = 
        vtkDoubleArray::SafeDownCast(
          qualityFilter->GetOutput()->GetFieldData()->GetArray("Mesh Tetrahedron Quality"));
      break;
    }
    case(3):
    {
      qualityField = 
        vtkDoubleArray::SafeDownCast(
          qualityFilter->GetOutput()->GetFieldData()->GetArray("Mesh Hexahedron Quality"));
      break;
    }
    default:
    {
      std::cout << "Invalid cell type. Only tri,quad,tet and hex supported." << std::endl;
      exit(1);
    }
  } 

  return qualityField;
}



void MeshQuality::cfmOptimize()
{
#ifdef HAVE_CFMSH

    // initialize openfoam
    //#include "setRootCase.H"
    int argc = 1;
    char** argv = new char*[2];
    argv[0] = new char[100];
    strcpy(argv[0], "NONE");
    Foam::argList args(argc, argv);
    //#include "createTime.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);
    //- 2d cartesian mesher cannot be run in parallel
    Foam::argList::noParallel();  


    //- load the mesh from disk
    polyMeshGen pmg(runTime);
    pmg.read();

    //- construct the smoother
    Foam::meshOptimizer mOpt(pmg);


    if( _cfmQPrms->_withConstraint )
    {
      std::string constrainedCellSet = _cfmQPrms->consCellSet;

      //- lock cells in constrainedCellSet
      mOpt.lockCellsInSubset(constrainedCellSet);

      //- find boundary faces which shall be locked
      Foam::labelLongList lockedBndFaces, selectedCells;

      const label sId = pmg.cellSubsetIndex(constrainedCellSet);
      pmg.cellsInSubset(sId, selectedCells);

      Foam::boolList activeCell(pmg.cells().size(), false);
      forAll(selectedCells, i)
          activeCell[selectedCells[i]] = true;
    }

    //- clear geometry information before volume smoothing
    pmg.clearAddressingData();

    //- perform optimisation using the laplace smoother and
    mOpt.optimizeMeshFV
    (
        _cfmQPrms->nLoops,
        _cfmQPrms->nLoops,
        _cfmQPrms->nIterations,
        _cfmQPrms->nSrfItr
    );

    //- perform optimisation of worst quality faces
    mOpt.optimizeMeshFVBestQuality(_cfmQPrms->nLoops, _cfmQPrms->qualThrsh);

    //- check the mesh again and untangl bad regions if any of them exist
    mOpt.untangleMeshFV(_cfmQPrms->nLoops, _cfmQPrms->nIterations, 
            _cfmQPrms->nSrfItr);
    Foam::Info << "Writing mesh" << endl;
    pmg.write();

    delete[] argv;
    Foam::Info << "Finished mesh quality enhancement task." << endl;
#else
    std::cerr << "Compile with CFMSH option enabled to use this feature.\n";
#endif
}

