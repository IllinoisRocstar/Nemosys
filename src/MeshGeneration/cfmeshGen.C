#include <iostream>
#include <string>
#include <cfmeshGen.H>
#include <cfmeshParams.H>
#include <boost/filesystem.hpp>

// vtk
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>
#include <vtkCellArray.h>

// openfoam headers
#include "fvCFD.H"
#include "fvMesh.H"
#include "vtkTopo.H"
#include "fileName.H"

// cfmesh headers
#include "cartesian2DMeshGenerator.H"
#include "triSurfaceDetectFeatureEdges.H"
#include "triSurfacePatchManipulator.H"
#include "triSurf.H"
#include "tetMeshGenerator.H"
#include "polyMeshGenModifier.H"
#include "meshOptimizer.H"
#include "cartesianMeshGenerator.H"
#include "voronoiMeshGenerator.H"


cfmeshGen::cfmeshGen() 
{
  // default meshing parameters
  _params = new cfmeshParams();
  defaults = true;

  // Initialization tasks
  initialize();
}

cfmeshGen::cfmeshGen(cfmeshParams* params):
    defaults(false), _params(params)
{
    // Initialization tasks
    initialize();
}

cfmeshGen::~cfmeshGen()
{
    if (defaults)
        delete _params;
}


void cfmeshGen::initialize()
{
    // surface feature edge treatment
    if (_params->srfEdge.has_value())
        if ( surfaceFeatureEdgeDetect() )
        {
            std::cerr << "A problem occured during edge detection step!\n";
            throw;
        }


    // create dictionaries needed
    createControlDict();
    createMshDict();
    createfvSchemesDict();
    createfvSolutionDict();

    // initialize openfoam
    //#include "setRootCase.H"
    int argc = 1;
    char** argv = new char*[2];
    argv[0] = new char[100];
    strcpy(argv[0], "NONE");
    _args = new Foam::argList(argc, argv);
    //#include "createTime.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    _runTime = new Foam::Time(Foam::Time::controlDictName, *_args);
    //- 2d cartesian mesher cannot be run in parallel
    Foam::argList::noParallel();  
}


void cfmeshGen::createControlDict()
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
deltaT	1;\n\n\
startTime	0;\n\n\
writeInterval	1;\n\n\
// ************************************************************************* //";
    contDict << contText;
    contDict.close();
}

void cfmeshGen::createfvSchemesDict()
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
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
// ************************************************************************* //";
    contDict << contText;
    contDict.close();
}

void cfmeshGen::createfvSolutionDict()
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
\n\
// ************************************************************************* //";
    contDict << contText;
    contDict.close();
}


void cfmeshGen::createMshDict()
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

    // creating mesh dictionary file
    std::ofstream contDict;
    contDict.open(std::string(dir_path)+"/meshDict");
    // header
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
    object    meshDict;\n\
}\n\n";

    contText = contText + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    
    if ((_params->maxCellSize)>0)
        contText = contText + "\nmaxCellSize " + std::to_string(_params->maxCellSize) + ";\n";
   
    if ((_params->minCellSize)>0)
        contText = contText + "\nminCellSize " + std::to_string(_params->minCellSize) + ";\n";

    if ((_params->bndryCellSize)>0)
        contText = contText + "\nboundaryCellSize " + 
        std::to_string(_params->bndryCellSize) + ";\n";

    if ((_params->bndryCellSizeRefThk)>0)
        contText = contText + "\nboundaryCellSizeRefinementThickness " + 
        std::to_string(_params->bndryCellSizeRefThk) + ";\n";

    
    if ((_params->keepCellIB) == true)
        contText = contText + "\nkeepCellsIntersectingBoundary 1;\n"; 
        //+ std::to_string(_params->keepCellIB) + ";\n";

    if ((_params->chkGluMsh) == true)
        contText = contText + "\ncheckForGluedMesh 1;\n"; 
        //+ std::to_string(_params->chkGluMsh) + ";\n";

    if ((_params->alwDiscDomains))
        contText = contText + "\nallowDisconnectedDomains 1;\n";

    contText = contText +"\n\nsurfaceFile	\"" + (_params->geomFilePath) +"\";\n";


    // boundary layer
    if (_params->boundaryLayers.has_value())
    {
        const auto &boundaryLayer = _params->boundaryLayers.value();
        contText = contText + "\nboundaryLayers\n{\n";
        contText = contText + "\tnLayers\t" + std::to_string(boundaryLayer.blNLyr)
            + ";\n";
        contText = contText + "\tthicknessRatio\t" + std::to_string(boundaryLayer.blThkRto) + ";\n";
        contText = contText + 
            ( (boundaryLayer.maxFrstLyrThk) > 0 ? "\n" : ("\tmaxFirstLayerThickness\t" +
              std::to_string(boundaryLayer.maxFrstLyrThk) + ";\n") );
        contText = contText + 
            ( (boundaryLayer.alwDiscont) ? ("\tallowDiscontinuity\t1;\n") : "\n" );
        
        // boundary layer patches
        if (!boundaryLayer.blPatches.empty())
        {
            contText = contText + "\tpatchBoundaryLayers\n\t{\n";
            for (auto pt=(boundaryLayer.blPatches).begin(); pt!=(boundaryLayer.blPatches).end(); pt++)
            {
                contText = contText + "\t\t\"" + (pt->patchName) + "\"\n\t\t{\n";
                if ( (pt->alwDiscont) == true)
                    contText = contText + "\t\t\tallowDiscontinuity\t1;\n";
                if ( (pt->maxFrstLyrThk) > 0)
                    contText = contText + "\t\t\tmaxFirstLayerThickness\t"
                        + std::to_string((pt->maxFrstLyrThk)) + ";\n";
                if ( (pt->blNLyr) > 0)
                    contText = contText + "\t\t\tnLayers\t"
                        + std::to_string((pt->blNLyr)) + ";\n";
                if ( (pt->blThkRto) > 0)
                    contText = contText + "\t\t\tthicknessRatio\t"
                        + std::to_string((pt->blThkRto)) + ";\n";
                contText = contText + "\t\t}\n";
            }
            contText = contText + "\n\t}\n";
        }

        contText = contText + "\n}\n\n";
    }

    // object refinements
    if (!_params->objRefLst.empty())
    {
        contText = contText + "objectRefinements\n{\n";
        for (auto ref=(_params->objRefLst).begin(); ref!=(_params->objRefLst).end(); ref++)
        {
            contText = contText + "\t" + (ref->name) +"\n\t{\n";
            for (auto prm=(ref->params).begin(); prm!=(ref->params).end(); prm++)
                contText = contText + "\t\t" + (prm->first) + "\t" + (prm->second) + ";\n";
            contText = contText + "\t}\n";
        }
        contText = contText+ "}\n";
    }

    // local refinement
    if (!_params->refPatches.empty())
    {
        contText = contText + "localRefinement\n{\n";
        for (auto pt=(_params->refPatches).begin(); pt!=(_params->refPatches).end(); pt++)
        {
            contText = contText + "\t\"" + (pt->patchName) + "\"\n\t{\n";
            if ( (pt->aditRefLvls) > 0)
                contText = contText + "\t\tadditionalRefinementLevels\t"
                    + std::to_string((pt->aditRefLvls)) + ";\n";
            if ( (pt->refThickness) > 0)
                contText = contText + "\t\trefinementThickness\t"
                    + std::to_string((pt->refThickness)) + ";\n";
            if ( (pt->cellSize) > 0)
                contText = contText + "\t\tcellSize\t"
                    + std::to_string((pt->cellSize)) + ";\n";
            contText = contText + "\t}\n";
        }
        contText = contText + "\n}\n";
    }

    // rename boundaries
    if (_params->renBndry.has_value())
    {
        const auto &renBndry = _params->renBndry.value();
        contText = contText + "renameBoundary\n{\n";
        contText = contText + "\tdefaultName\t" + (renBndry).defName +";\n";
        contText = contText + "\tdefaultType\t" + (renBndry).defType +";\n";
        contText = contText + "\tnewPatchNames\n\t{\n";
        for (auto pt=(renBndry).newPatches.begin();
                pt!=(renBndry).newPatches.end(); pt++)
        {
            contText = contText + "\t\t\"" + (pt->name) + "\"\n\t\t{\n";
            contText = contText + "\t\t\tnewName\t" + pt->newName + ";\n";
            contText = contText + "\t\t\ttype\t" + pt->newType + ";\n";
            contText = contText + "\t\t}\n";
        }
        contText = contText + "\t\n}\n";
        contText = contText + "\n}\n";
    }

    contText = contText +
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
    contDict << contText;
    contDict.close();
}


int cfmeshGen::createMeshFromSTL(const char* fname)
{
    // mesh generation and I/O
    Foam::Info << "Generating mesh with cfMesh engine" << Foam::endl;
    if (_params->generator == "cartesian2D")
    {
        Foam::cartesian2DMeshGenerator cmg(*_runTime);
        std::cout << "ExecutionTime = " << _runTime->elapsedCpuTime() << " s\n"
            << "ClockTime = " << _runTime->elapsedClockTime() << " s" << std::endl;
        cmg.writeMesh();
    }
    else if (_params->generator == "tetMesh")
    {
        Foam::tetMeshGenerator tmg(*_runTime);
        std::cout << "ExecutionTime = " << _runTime->elapsedCpuTime() << " s\n"
            << "ClockTime = " << _runTime->elapsedClockTime() << " s" << std::endl;
        tmg.writeMesh();

        // post-processing steps
        if (_params->improveMeshQuality.has_value())
            improveMeshQuality();
    }
    else if (_params->generator == "cartesian3D")
    {
        Foam::cartesianMeshGenerator cmg(*_runTime);
        std::cout << "ExecutionTime = " << _runTime->elapsedCpuTime() << " s\n"
            << "ClockTime = " << _runTime->elapsedClockTime() << " s" << std::endl;
        cmg.writeMesh();
    }
    else if (_params->generator == "polyMesh")
    {
        Foam::voronoiMeshGenerator pmg(*_runTime);
        std::cout << "ExecutionTime = " << _runTime->elapsedCpuTime() << " s\n"
            << "ClockTime = " << _runTime->elapsedClockTime() << " s" << std::endl;
        pmg.writeMesh();
    } 
    else
    {
        std::cerr 
            << (_params->generator)
            << " is not a supported mesh generator.\n";
        throw;
    }

    // loading the mesh file
    readFoamMesh();
}

void cfmeshGen::readFoamMesh()
{
    // reading mesh database and converting
    Foam::word regionName;
    if (_args->optionReadIfPresent("region", regionName))
    {
        Foam::Info
            << "Create mesh " << regionName << " for time = "
            << _runTime->timeName() << Foam::nl << Foam::endl;
    }
    else
    {
        regionName = Foam::fvMesh::defaultRegion;
        Foam::Info
            << "Create mesh for time = "
            << _runTime->timeName() << Foam::nl << Foam::endl;
    }
    _fmesh = new Foam::fvMesh 
    (
        Foam::IOobject
        (
            regionName,
            _runTime->timeName(),
            *_runTime,
            Foam::IOobject::MUST_READ
        )
    );  

    genMshDB();
}

void cfmeshGen::genMshDB()
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
            dataSet_tmp->InsertNextCell(VTK_POLYHEDRON, pntIdsVec.size(), &pntIdsVec[0],
                    nFace, faces->GetPointer());
        }
    }
    dataSet = dataSet_tmp;
}

void cfmeshGen::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                              const int cellType, 
                              std::vector<int>& vrtIds)
{
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(vrtIds.size());
    for (auto pit = vrtIds.begin(); pit!= vrtIds.end(); pit++)
        vtkCellIds->SetId(pit-vrtIds.begin(), *pit);
    dataSet->InsertNextCell(cellType,vtkCellIds);
}
                             
int cfmeshGen::surfaceFeatureEdgeDetect()
{
    std::cout << "Performing surface feature edge detection.\n";
    std::string of = "./"+caseName+"_feature.ftr";
    Foam::fileName inFileName(_params->geomFilePath);
    fileName outFileName(of);

    if (outFileName == inFileName)
    {
        std::cerr << "Output file " << outFileName
            << " would overwrite the input file.\n";
        throw;
    }

    double tol = _params->srfEdge.value().srfEdgAng;
    std::cout << "Using " << tol <<" deg angle\n";

    Foam::triSurf originalSurface(inFileName);

    Foam::triSurfaceDetectFeatureEdges edgeDetector(originalSurface, tol);
    edgeDetector.detectFeatureEdges();

    if( outFileName.ext() == "fms" || outFileName.ext() == "FMS" )
    {
        std::cout << "Writing : " << outFileName << std::endl;
        originalSurface.writeSurface(outFileName);
    }
    else
    {
        Foam::triSurfacePatchManipulator manipulator(originalSurface);
        const triSurf* newSurfPtr = manipulator.surfaceWithPatches();

        std::cout << "Writing : " << outFileName << std::endl;
        newSurfPtr->writeSurface(outFileName);

        delete newSurfPtr;
    }

    // change cad input file 
    _params->geomFilePath = of;

    return 0;
}

int cfmeshGen::improveMeshQuality()
{
    std::cout << "Performing mesh quality improvements.\n";

    //- load the mesh from disk
    Foam::polyMeshGen pmg(*_runTime);
    pmg.read();

    //- construct the smoother
    Foam::meshOptimizer mOpt(pmg);

    const auto &meshQual = _params->improveMeshQuality.value();

    if( (meshQual.qltConCelSet) != "none" )
    {
        //- lock cells in constrainedCellSet
        mOpt.lockCellsInSubset((meshQual.qltConCelSet));

        //- find boundary faces which shall be locked
        Foam::labelLongList lockedBndFaces, selectedCells;

        const Foam::label sId = pmg.cellSubsetIndex((meshQual.qltConCelSet));
        pmg.cellsInSubset(sId, selectedCells);

        Foam::boolList activeCell(pmg.cells().size(), false);
        for (int iCl=0; iCl< selectedCells.size(); iCl++)
            activeCell[selectedCells[iCl]] = true;
    }

    //- clear geometry information before volume smoothing
    pmg.clearAddressingData();

    //- perform optimisation using the laplace smoother and
    mOpt.optimizeMeshFV
    (
        (meshQual.qltNLop),
        (meshQual.qltNLop),
        (meshQual.qltNItr),
        (meshQual.qltNSrfItr)
    );

    //- perform optimisation of worst quality faces
    mOpt.optimizeMeshFVBestQuality((meshQual.qltNLop), (meshQual.qltQltThr));

    //- check the mesh again and untangl bad regions if any of them exist
    mOpt.untangleMeshFV((meshQual.qltNLop), (meshQual.qltNItr), (meshQual.qltNSrfItr));

    std::cout << "Finished optimization cycle\n";
    pmg.write();

    return 0;    
}

