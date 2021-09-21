#include "AuxiliaryFunctions.H"
#include <vtkAppendFilter.h>
#include <vtkCell.h>
#include <vtkCell3D.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkUnstructuredGrid.h>

#include <boost/filesystem.hpp>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "MeshManipulationFoam/MeshManipulationFoam.H"

// New
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

// openfoam headers
#include <fileName.H>
#include <fvMesh.H>
#include <foamVtkVtuAdaptor.H>
#include <getDicts.H>
#include <surfaceLambdaMuSmooth.H>  // SurfLambdaMuSmooth
#include <splitMeshRegions.H>       // splitMeshByRegions
#include <mergeMeshes.H>            // mergeMeshes
#include <createPatch.H>            // createPatch
#include <surfaceSplitByTopology.H> // surfaceSplitByTopology
#include <foamToSurface.H>          // foamToSurface

// third party
#include <ANN/ANN.h>

// TODO
// 1. Two of these methods (mergeMesh and CreatePatch) are little bit pack
//    mesh specific but can be modified to accept very broad user
//    arguments and perform mesh manipulations.

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

MeshManipulationFoam::~MeshManipulationFoam() {
  // class destructor
}

/*
  SurfaceLambdaMuSmooth takes input surface file and user defined values of
  lambda, mu, and surface smoothing iterations to perform laplacian smoothing
  on surface of geometry. This utility accepts following extensions as input:
  .ofs, .obj, .inp, .stl, .tri, .off, .stlb, .nas, .bdf, .gts, .vtk, and .ac.
  It will write output file at userdefined path
*/
void MeshManipulationFoam::surfLambdaMuSmooth() {
  const auto &surfLMSmoothParams = _mshMnipPrms->surfLMSmoothParams;
  const Foam::fileName surfFileName = (surfLMSmoothParams.slmssurfaceFile);
  const Foam::fileName outFileName = (surfLMSmoothParams.slmsoutputFile);
  const Foam::scalar lambda = (surfLMSmoothParams.lambda_);
  const Foam::scalar mu = (surfLMSmoothParams.mu);
  const Foam::label iters = (surfLMSmoothParams.slmsIterations);
  const bool addFtrFl = (surfLMSmoothParams.addFeatureFile);
  auto slmsObj = surfaceLambdaMuSmooth();
  slmsObj.execute(surfFileName, outFileName, lambda, mu, iters, addFtrFl);
}

/*
  SplitMeshRegions utility splits mesh into multiple regions. It takes cell-face
  -cell walk through the mesh, separates different cellZones into domains, and
  writes them into constant directory. Currently it supports splitting of mesh
  using cellZones. This function returns an integer to be used in mergeMeshes
  utility. SplitMeshRegions utility assigns random numbering to split domains
  from 0 to number of domains and skips a number automatically. This number is
  returned and used in mergeMeshes function so that it knows which directory is
  missing from constant folder.
*/
std::pair<std::vector<int>, std::vector<std::string>>
MeshManipulationFoam::splitMshRegions() {
  const auto &spltMshRegions = _mshMnipPrms->splitMeshRegParams;
  const bool useCellZones = (spltMshRegions.cellZones);
  auto smrObj = splitMeshRegions();
  return smrObj.execute(useCellZones);
}

/*
MergeMeshes utility takes master region and slave region names/paths from user
and merges the slave regions to master region one by one.
*/
void MeshManipulationFoam::mergeMesh(int dirStat, int nDomains) {
  auto mmObj = mergeMeshes();
  const auto &mergeMeshesParams = _mshMnipPrms->mergeMeshesParams;

  std::vector<std::string> addCases;

  // Collecting all domain names for automatic merging of all mesh regions.
  // Directory number obtained from splitMeshRegions is used here. Also
  // number of packs are passed through this function for use in loop.

  if (nDomains == 2) {
    addCases.push_back(mergeMeshesParams.addCase);
  } else {
    for (int i = 2; i < (nDomains + 1); i++) {
      if (i == dirStat) i++;
      addCases.push_back("domain" + (std::to_string(i)));
    }
    addCases.push_back(mergeMeshesParams.addCase);
  }
  mmObj.execute(mergeMeshesParams.masterCasePath, mergeMeshesParams.addCasePath,
                addCases, nDomains, dirStat);
}

/*
CreatePatch utility merges multiple patches of certain mesh in a single patch.
It reads createPatchDict from system/mesh_reg_name/createPatchDict. Currently
the utility implementation is wrapped with for loop to execute it twice for
Pack Mesh Generation service.
*/
void MeshManipulationFoam::createPtch(int dirStat) {
  const auto &createPatchParams = _mshMnipPrms->createPatchParams;
  auto cpObj = createPatch();
  createPatchDict(dirStat, true);
  cpObj.execute(dirStat, createPatchParams.pathSurrounding, cpDict_);
}

/*
This utility converts OpenFoam mesh into STL file. Filename is user-input.
It can take paths for output file location and name.
*/
void MeshManipulationFoam::foamToSurf() {
  const auto &foamToSurfaceParams = _mshMnipPrms->foamToSurfParams;
  auto ftsObj = foamToSurface();
  ftsObj.execute(foamToSurfaceParams.outSurfName);
}

/*
surfaceSplitByTopology is a surface file manipulation utility which aims at
splitting multiple disconnected regions in geometry into separate surfaces and
outputs a single surface file containing all regions as well as separate STL
files for all regions. User inputs are input file name and output file name.
*/
int MeshManipulationFoam::surfSpltByTopology() {
  const auto &surfSplitParams = _mshMnipPrms->surfSplitParams;
  Foam::fileName surfFileName(surfSplitParams.surfFile);
  Foam::Info << "Reading surface from " << surfFileName << Foam::endl;
  Foam::fileName outFileName(surfSplitParams.outSurfFile);
  auto ssbtObj = surfaceSplitByTopology();
  return ssbtObj.execute(surfFileName, outFileName);
}

/*
CreatePatchDict function creates dictionary file for createPatch utility and
puts them in defined path (i.e system/domainX). Currently this function has
customizations for Pack Mesh Generation service. More general methods could be
added in future.
*/
void MeshManipulationFoam::createPatchDict(const int &dirStat,
                                           const bool &write) {
  const auto &createPatchParams = _mshMnipPrms->createPatchParams;

  // creating a base system directory
  std::string dir_path = "./system/" + createPatchParams.pathSurrounding;

  if (write) {
    boost::filesystem::path dir(dir_path);
    try {
      boost::filesystem::create_directory(dir);
    } catch (boost::filesystem::filesystem_error &e) {
      std::cerr << "Problem in creating system directory for the snappyHexMesh"
                << "\n";
      std::cerr << e.what() << std::endl;
      throw;
    }
  }

  // header
  std::string contText =
      "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: snappyHexMesh interface               |\n\
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
    object    createPatchDict;\n\
}\n\n";

  contText =
      contText +
      "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

  contText = contText + "\n\npointSync true;\n";
  contText = contText + "\n\npatches\n";
  contText = contText + "(\n";
  contText = contText + "\t{\n";
  contText =
      contText + "\t\tname " + (createPatchParams.surroundingName) + ";\n";
  contText = contText + "\n\t\tpatchInfo\n";
  contText = contText + "\t\t{\n";
  contText =
      contText + "\t\t\ttype " + (createPatchParams.srrndngPatchType) + ";\n";
  contText = contText + "\t\t}\n";
  contText = contText + "\n\t\tconstructFrom patches;\n";
  contText = contText + "\n\t\tpatches (\"domain0_to_domain.*\");";
  contText = contText + "\n\t}\n";
  contText = contText + ");\n";

  contText = contText +
             "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
             "*//\n\n";

  // creating mesh dictionary file
  std::ofstream contDict;
  if (write) {
    contDict.open(std::string(dir_path) + "/createPatchDict");
    contDict << contText;
    contDict.close();
  }

  // Keep this dictionary in memory
  Foam::dictionary tmptmpDuc =
      new Foam::dictionary(Foam::IStringStream(contText)(), true);
  cpDict_ = std::unique_ptr<Foam::dictionary>(
      new Foam::dictionary("createPatchDict"));
  cpDict_->merge(tmptmpDuc);

  if (dirStat == 1) {
    const char dir_path2[] = "./system/domain2";

    if (write) {
      boost::filesystem::path dir2(dir_path2);
      try {
        boost::filesystem::create_directory(dir2);
      } catch (boost::filesystem::filesystem_error &e) {
        std::cerr << "Problem in creating system directory for createPatch"
                  << "\n";
        std::cerr << e.what() << std::endl;
        throw;
      }
    }

    // header
    std::string contText2 =
        "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: snappyHexMesh interface               |\n\
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
    object    createPatchDict;\n\
}\n\n";

    contText2 = contText2 +
                "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
                "* * //\n\n";

    contText2 = contText2 + "\n\npointSync true;\n";
    contText2 = contText2 + "\n\npatches\n";
    contText2 = contText2 + "(\n";
    contText2 = contText2 + "\t{\n";
    contText2 = contText2 + "\t\tname " + (createPatchParams.packsName) + ";\n";
    contText2 = contText2 + "\n\t\tpatchInfo\n";
    contText2 = contText2 + "\t\t{\n";
    contText2 =
        contText2 + "\t\t\ttype " + (createPatchParams.packsPatchType) + ";\n";
    contText2 = contText2 + "\t\t}\n";
    contText2 = contText2 + "\n\t\tconstructFrom patches;\n";
    contText2 = contText2 + "\n\t\tpatches (\"domain.*_to_domain0\");";
    contText2 = contText2 + "\n\t}\n";
    contText2 = contText2 + ");\n";

    contText2 = contText2 +
                "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
                "* * //\n\n";

    // creating mesh dictionary file
    std::ofstream contDict2;
    if (write) {
      contDict2.open(std::string(dir_path2) + "/createPatchDict");
      contDict2 << contText2;
      contDict2.close();
    }

    // Keep this dictionary in memory
    cpDict_.reset();
    Foam::dictionary tmptmpDuc2 =
        new Foam::dictionary(Foam::IStringStream(contText2)(), true);
    cpDict_ = std::unique_ptr<Foam::dictionary>(
        new Foam::dictionary("createPatchDict"));
    cpDict_->merge(tmptmpDuc2);
  } else {
    const char dir_path2[] = "./system/domain1";
    if (write) {
      boost::filesystem::path dir2(dir_path2);
      try {
        boost::filesystem::create_directory(dir2);
      } catch (boost::filesystem::filesystem_error &e) {
        std::cerr << "Problem in creating system directory for createPatch"
                  << "\n";
        std::cerr << e.what() << std::endl;
        throw;
      }
    }

    // header
    std::string contText2 =
        "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: snappyHexMesh interface               |\n\
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
    object    createPatchDict;\n\
}\n\n";

    contText2 = contText2 +
                "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
                "* * //\n\n";

    contText2 = contText2 + "\n\npointSync true;\n";
    contText2 = contText2 + "\n\npatches\n";
    contText2 = contText2 + "(\n";
    contText2 = contText2 + "\t{\n";
    contText2 = contText2 + "\t\tname " + (createPatchParams.packsName) + ";\n";
    contText2 = contText2 + "\n\t\tpatchInfo\n";
    contText2 = contText2 + "\t\t{\n";
    contText2 =
        contText2 + "\t\t\ttype " + (createPatchParams.packsPatchType) + ";\n";
    contText2 = contText2 + "\t\t}\n";
    contText2 = contText2 + "\n\t\tconstructFrom patches;\n";
    contText2 = contText2 + "\n\t\tpatches (\"domain.*_to_domain0\");";
    contText2 = contText2 + "\n\t}\n";
    contText2 = contText2 + ");\n";

    contText2 = contText2 +
                "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
                "* * //\n\n";

    // creating mesh dictionary file
    std::ofstream contDict2;
    if (write) {
      contDict2.open(std::string(dir_path2) + "/createPatchDict");
      contDict2 << contText2;
      contDict2.close();
    }

    // Keep this dictionary in memory
    Foam::dictionary tmptmpDuc2 =
        new Foam::dictionary(Foam::IStringStream(contText2)(), true);
    cpDict_.reset();
    cpDict_ = std::unique_ptr<Foam::dictionary>(
        new Foam::dictionary("createPatchDict"));
    cpDict_->merge(tmptmpDuc2);
  }
}

// Adds connectivity information at shared interfaces.
void MeshManipulationFoam::addCohesiveElements(double tol,
                                               const std::string &outName) {
  bool writeDicts = true;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);

  Foam::Time runTime(controlDict_.get(), ".", ".");

  auto _fmeshSurr = std::unique_ptr<NEM::MSH::foamGeoMesh>(
      NEM::MSH::foamGeoMesh::Read("domain0.foam"));

  auto _fmeshPacks = std::unique_ptr<NEM::MSH::foamGeoMesh>(
      NEM::MSH::foamGeoMesh::Read("domain1.foam"));

  // declare vtk datasets
  vtkSmartPointer<vtkUnstructuredGrid> dataSetSurr =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkUnstructuredGrid> surrTmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto objVfoamSurr = Foam::vtk::vtuAdaptor();
  surrTmp = objVfoamSurr.internal(_fmeshSurr->getFoamMesh());

  int impPts = (int)(surrTmp->GetNumberOfPoints());

  vtkSmartPointer<vtkUnstructuredGrid> pckTmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto objVfoamPck = Foam::vtk::vtuAdaptor();
  pckTmp = objVfoamPck.internal(_fmeshPacks->getFoamMesh());

  // Merge two foam meshes
  vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
  appendFilter->AddInputData(surrTmp);
  appendFilter->AddInputData(pckTmp);
  appendFilter->MergePointsOn();
  appendFilter->Update();
  dataSetSurr = appendFilter->GetOutput();

  int numPoints = (int)(dataSetSurr->GetNumberOfPoints());

  int nDim = 3;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(numPoints, nDim);
  for (int iPnt = 0; iPnt < numPoints; iPnt++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iPnt, &getPt[0]);
    pntCrd[iPnt][0] = getPt[0];
    pntCrd[iPnt][1] = getPt[1];
    pntCrd[iPnt][2] = getPt[2];
  }

  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, numPoints, nDim);

  // Finding duplicate nodes
  double rad = 1e-05;
  std::map<int, int> dupNdeMap;
  ANNpoint qryPnt;  // Query point.
  int nNib = 2;  // Number of neighbours to return (including query pt. itself).
  ANNidxArray nnIdx = new ANNidx[nNib];    // Nearest neighbour ID array.
  ANNdistArray dists = new ANNdist[nNib];  // Neighbour distance array.
  qryPnt = annAllocPt(nDim);               // Initializing query point

  for (int iNde = 0; iNde < impPts; iNde++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iNde, &getPt[0]);
    qryPnt[0] = getPt[0];
    qryPnt[1] = getPt[1];
    qryPnt[2] = getPt[2];
    kdTree->annkFRSearch(qryPnt, rad, nNib, nnIdx, dists, 0);
    if (dists[1] <= tol) {
      if ((nnIdx[1]) >= impPts)
        dupNdeMap[nnIdx[0]] = nnIdx[1];  // Format Surrounding::Pack
      else
        dupNdeMap[nnIdx[1]] = nnIdx[0];  // Format Surrounding::Pack
    }
  }

  int numDupPnts = (int)dupNdeMap.size();
  std::cout << "Found " << numDupPnts << " duplicate nodes.\n";

  // Getting cells in packs
  std::map<int, int>::iterator it = dupNdeMap.begin();
  std::vector<int> cellID;
  for (int i = 0; i < (int)(dupNdeMap.size()); i++) {
    int nCells, cellNum;
    vtkIdList *cells = vtkIdList::New();

    // get cells the point belongs to
    dataSetSurr->GetPointCells(it->second, cells);
    nCells = (int)cells->GetNumberOfIds();

    for (cellNum = 0; cellNum < nCells; cellNum++) {
      cellID.push_back((int)(cells->GetId(cellNum)));  // get cell id from list
    }
    it++;
  }

  sort(cellID.begin(), cellID.end());
  cellID.erase(unique(cellID.begin(), cellID.end()), cellID.end());

  int know = (int)cellID.size();

  // std::cout << "Total Cells Found Are " << know << std::endl;

  // Determines which cells are useful
  std::vector<int> newCellIds;
  std::vector<int> debugCellIds;

  for (int i = 0; i < know; i++) {
    // Getting all cell defining points
    std::vector<int> ptIdzz(8);
    int cntr = 0;
    vtkCell *cell;
    cell = dataSetSurr->GetCell(cellID[i]);

    vtkIdList *pts = cell->GetPointIds();

    for (int j = 0; j < 8; j++) { ptIdzz[j] = (int)(pts->GetId(j)); }

    // Checking if any cell contains any 4 or more of duplicate nodes.
    std::map<int, int>::iterator it2 = dupNdeMap.begin();
    while (it2 != dupNdeMap.end()) {
      for (int k = 0; k < 8; k++) {
        if ((it2->second) == ptIdzz[k]) {
          cntr++;
        } else {
          // Nothing
        }
      }
      it2++;
    }

    if (cntr >= 4) {
      newCellIds.push_back(cellID[i]);
    } else {
      debugCellIds.push_back(cellID[i]);
    }
  }

  int size2 = (int)newCellIds.size();

  // Getting useful faces from cells
  std::vector<int> globalPtIds;
  std::vector<int> surroundingArray;
  std::vector<int> packArray;
  for (int i = 0; i < size2; i++) {
    vtkCell *cell;
    cell = dataSetSurr->GetCell(newCellIds[i]);
    vtkIdList *pts = cell->GetPointIds();

    for (int j = 0; j < 6; j++) {
      int isFour = 0;
      int *ptFaces = nullptr;
      vtkCell3D *cell3d = static_cast<vtkCell3D *>(cell);
      cell3d->GetFacePoints(j, ptFaces);
      std::vector<int> keysDupMap = std::vector<int>(4);
      std::map<int, int>::iterator it3 = dupNdeMap.begin();
      while (it3 != dupNdeMap.end()) {
        for (int h = 0; h < 4; h++) {
          if ((it3->second) == (pts->GetId(ptFaces[h]))) {
            isFour++;
            keysDupMap[h] = it3->first;
          } else {
            // Nothing
          }
        }

        it3++;
      }

      if (isFour == 4) {
        for (int k = 0; k < 4; k++) { globalPtIds.push_back(keysDupMap[k]); }
      }
    }
  }

  std::map<int, int>::iterator it5;

  // Creating node map with sequence.
  std::unordered_multimap<int, int> cohesiveMap;

  for (int i = 0; i < (int)globalPtIds.size(); i++) {
    it5 = dupNdeMap.find(globalPtIds[i]);

    if (it5 == dupNdeMap.end()) {
      // Nothing
    } else {
      surroundingArray.push_back(it5->first);
      packArray.push_back(it5->second);
    }
  }

  // std::cout << "Size = " << packArray.size() << std::endl;
  int newCells = (int)(packArray.size() / 4);
  int startNum = (int)(dataSetSurr->GetNumberOfCells());

  // Finally, creating VTK cells
  int it7 = 0;
  std::vector<int> pntCohesiveIds;
  int nCelCohesivePnts;
  for (int icl = 0; icl < newCells; icl++) {
    nCelCohesivePnts = 4;
    pntCohesiveIds.resize((nCelCohesivePnts * 2), -1);
    for (int ip = 0; ip < nCelCohesivePnts; ip++) {
      pntCohesiveIds[ip] = packArray[it7];
      pntCohesiveIds[ip + 4] = surroundingArray[it7];
      it7++;
    }
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(pntCohesiveIds.size());
    for (auto pit = pntCohesiveIds.begin(); pit != pntCohesiveIds.end(); pit++)
      vtkCellIds->SetId(pit - pntCohesiveIds.begin(), *pit);
    dataSetSurr->InsertNextCell(12, vtkCellIds);
  }

  int endNum = (int)(dataSetSurr->GetNumberOfCells());

  // ****************************** //
  // Zero volume cells visualization method for debugging purpose.
  std::vector<int> cohCellIds;
  for (int i = startNum; i < endNum; i++) { cohCellIds.push_back(i); }

  int size3 = (int)cohCellIds.size();
  // Takes cell Ids as input and writes VTK mesh
  vtkSmartPointer<vtkUnstructuredGrid> visDataSet =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  std::vector<int> ptIdzz2;
  for (int i = 0; i < size3; i++) {
    vtkCell *cell;

    cell = dataSetSurr->GetCell(cohCellIds[i]);

    vtkIdList *pts = cell->GetPointIds();

    for (int j = 0; j < 8; j++) { ptIdzz2.push_back((int)(pts->GetId(j))); }
  }

  int lvar = 0;

  vtkSmartPointer<vtkPoints> pointsViz = vtkSmartPointer<vtkPoints>::New();

  for (int i = 0; i < (size3 * 8); i++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(ptIdzz2[i], &getPt[0]);
    pointsViz->InsertNextPoint(getPt[0], getPt[1], getPt[2]);
  }

  visDataSet->SetPoints(pointsViz);

  for (int i = 0; i < size3; i++) {
    std::vector<int> ptnewIdz;
    for (int j = 0; j < 8; j++) {
      ptnewIdz.push_back(lvar);
      lvar++;
    }
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(ptnewIdz.size());
    for (auto pit = ptnewIdz.begin(); pit != ptnewIdz.end(); pit++)
      vtkCellIds->SetId(pit - ptnewIdz.begin(), *pit);
    visDataSet->InsertNextCell(12, vtkCellIds);
  }

  auto vm = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
      new NEM::MSH::vtkGeoMesh(visDataSet));
  vm->write("CohesiveElements.vtu");

  // ********************************************* //
  auto vm2 = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
      new NEM::MSH::vtkGeoMesh(dataSetSurr));
  vm2->write(outName);

  // Convert final mesh to tetrahedral
  bool tetra = false;
  std::string ofnameTet = "Tetra" + outName;

  if (tetra) {
    vtkSmartPointer<vtkDataSetTriangleFilter> triFilter =
        vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triFilter->SetInputData(dataSetSurr);
    triFilter->Update();

    auto newData_tet = triFilter->GetOutput();

    auto vm3 = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
        new NEM::MSH::vtkGeoMesh(newData_tet));
    vm3->write(ofnameTet);
  }
  // **************************************************************************
}

void MeshManipulationFoam::addArtificialThicknessElements(
    double &tol, const std::string &outName, double &thickness) {
  // Initializing FOAM
  bool writeDicts = true;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);

  Foam::Time runTime(controlDict_.get(), ".", ".");

  auto _fmeshSurr = std::unique_ptr<NEM::MSH::foamGeoMesh>(
      NEM::MSH::foamGeoMesh::Read("domain0.foam"));

  auto _fmeshPacks = std::unique_ptr<NEM::MSH::foamGeoMesh>(
      NEM::MSH::foamGeoMesh::Read("domain1.foam"));

  // declare vtk datasets
  vtkSmartPointer<vtkUnstructuredGrid> dataSetSurr =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkUnstructuredGrid> surrTmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto objVfoamSurr = Foam::vtk::vtuAdaptor();
  surrTmp = objVfoamSurr.internal(_fmeshSurr->getFoamMesh());

  int impPts = (surrTmp->GetNumberOfPoints());

  vtkSmartPointer<vtkUnstructuredGrid> pckTmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto objVfoamPck = Foam::vtk::vtuAdaptor();
  pckTmp = objVfoamPck.internal(_fmeshPacks->getFoamMesh());

  // Merge two foam meshes
  vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
  appendFilter->AddInputData(surrTmp);
  appendFilter->AddInputData(surrTmp);
  appendFilter->MergePointsOn();
  appendFilter->Update();
  dataSetSurr = appendFilter->GetOutput();

  int numPoints = (int)(dataSetSurr->GetNumberOfPoints());

  int nDim = 3;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(numPoints, nDim);
  for (int iPnt = 0; iPnt < numPoints; iPnt++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iPnt, &getPt[0]);
    pntCrd[iPnt][0] = getPt[0];
    pntCrd[iPnt][1] = getPt[1];
    pntCrd[iPnt][2] = getPt[2];
  }

  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, numPoints, nDim);

  // Finding duplicate nodes
  double rad = 1e-05;
  std::map<int, int> dupNdeMap;
  ANNpoint qryPnt;  // Query point.
  int nNib = 2;  // Number of neighbours to return (including query pt. itself).
  ANNidxArray nnIdx = new ANNidx[nNib];    // Nearest neighbour ID array.
  ANNdistArray dists = new ANNdist[nNib];  // Neighbour distance array.
  qryPnt = annAllocPt(nDim);               // Initializing query point

  for (int iNde = 0; iNde < impPts; iNde++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iNde, &getPt[0]);
    qryPnt[0] = getPt[0];
    qryPnt[1] = getPt[1];
    qryPnt[2] = getPt[2];
    kdTree->annkFRSearch(qryPnt, rad, nNib, nnIdx, dists, 0);
    if (dists[1] <= tol) {
      if ((nnIdx[1]) >= impPts) {
        dupNdeMap[nnIdx[0]] = nnIdx[1];  // Format Surrounding::Pack
      } else {
        dupNdeMap[nnIdx[1]] = nnIdx[0];  // Format Surrounding::Pack
      }
    }
  }

  int numDupPnts = (int)dupNdeMap.size();
  std::cout << "Found " << numDupPnts << " duplicate nodes.\n";

  // Getting cells in packs
  std::map<int, int>::iterator it = dupNdeMap.begin();
  std::vector<int> cellID;
  for (int i = 0; i < (int)(dupNdeMap.size()); i++) {
    int nCells, cellNum;
    vtkIdList *cells = vtkIdList::New();

    // get cells the point belongs to
    dataSetSurr->GetPointCells(it->second, cells);
    nCells = (int)(cells->GetNumberOfIds());

    for (cellNum = 0; cellNum < nCells; cellNum++) {
      cellID.push_back(
          (int)(cells->GetId(cellNum)));  // get cell id from the list
    }
    it++;
  }

  sort(cellID.begin(), cellID.end());
  cellID.erase(unique(cellID.begin(), cellID.end()), cellID.end());

  int know = cellID.size();

  // Determines which cells are useful
  std::vector<int> newCellIds;
  std::vector<int> debugCellIds;

  for (int i = 0; i < know; i++) {
    // Getting all cell defining points
    std::vector<int> ptIdzz(8);
    int cntr = 0;
    vtkCell *cell;
    cell = dataSetSurr->GetCell(cellID[i]);

    vtkIdList *pts = cell->GetPointIds();

    for (int j = 0; j < 8; j++) { ptIdzz[j] = pts->GetId(j); }

    // Checking if any cell contains any 4 or more of duplicate nodes.
    std::map<int, int>::iterator it2 = dupNdeMap.begin();
    while (it2 != dupNdeMap.end()) {
      for (int k = 0; k < 8; k++) {
        if ((it2->second) == ptIdzz[k]) {
          cntr++;
        } else {
          // Nothing
        }
      }
      it2++;
    }

    if (cntr >= 4) {
      newCellIds.push_back(cellID[i]);
    } else {
      debugCellIds.push_back(cellID[i]);
    }
  }

  int size2 = (int)(newCellIds.size());

  // Getting useful faces from cells
  std::vector<int> globalPtIds;
  std::vector<int> surroundingArray;
  std::vector<int> packArray;
  for (int i = 0; i < size2; i++) {
    vtkCell *cell;
    cell = dataSetSurr->GetCell(newCellIds[i]);
    vtkIdList *pts = cell->GetPointIds();

    for (int j = 0; j < 6; j++) {
      int isFour = 0;
      int *ptFaces = nullptr;
      vtkCell3D *cell3d = static_cast<vtkCell3D *>(cell);
      cell3d->GetFacePoints(j, ptFaces);
      std::vector<int> keysDupMap = std::vector<int>(4);
      std::map<int, int>::iterator it3 = dupNdeMap.begin();
      while (it3 != dupNdeMap.end()) {
        for (int h = 0; h < 4; h++) {
          if ((it3->second) == (pts->GetId(ptFaces[h]))) {
            isFour++;
            keysDupMap[h] = it3->first;
          } else {
            // Nothing
          }
        }

        it3++;
      }

      if (isFour == 4) {
        for (int k = 0; k < 4; k++) { globalPtIds.push_back(keysDupMap[k]); }
      }
    }
  }

  std::map<int, int>::iterator it5;

  // Creating node map with sequence.
  std::unordered_multimap<int, int> cohesiveMap;

  for (int i = 0; i < (int)globalPtIds.size(); i++) {
    it5 = dupNdeMap.find(globalPtIds[i]);

    if (it5 == dupNdeMap.end()) {
      // Nothing
    } else {
      surroundingArray.push_back(it5->first);
      packArray.push_back(it5->second);
    }
  }

  std::cout << "End of previous method" << std::endl;

  // Converting vtkUnstructuredData to vtkPolyData
  vtkSmartPointer<vtkGeometryFilter> geometryFilter =
      vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputData(dataSetSurr);
  geometryFilter->Update();

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata = geometryFilter->GetOutput();

  vtkSmartPointer<vtkPolyDataNormals> polyDataNormal =
      vtkSmartPointer<vtkPolyDataNormals>::New();

  polyDataNormal->SetInputData(polydata);

  // Setting points normals calculation to ON with other options
  polyDataNormal->ComputePointNormalsOn();
  polyDataNormal->ComputeCellNormalsOff();
  // polyDataNormal->SetNonManifoldTraversal(1);
  // polyDataNormal->SetAutoOrientNormals(0);
  polyDataNormal->SetSplitting(0);
  // polyDataNormal->SetFeatureAngle(0.1);
  polyDataNormal->SetConsistency(1);
  // polyDataNormal->SetFlipNormals(0);
  polyDataNormal->Update();

  // Extracting point normals from polyDataNormal
  polydata = polyDataNormal->GetOutput();

  // Count points
  vtkIdType numPointsNew = polydata->GetNumberOfPoints();
  std::cout << "There are " << numPointsNew << " points." << std::endl;

  vtkDataArray *normalsGeneric = polydata->GetPointData()->GetNormals();

  // Number of normals should match number of points
  std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
            << " normals in normalsGeneric" << std::endl;

  // Changing point coordinates for adding artificial thickness;

  // Surrounding cells
  for (int i = 0; i < (int)surroundingArray.size(); i++) {
    double normalOfPt[3];
    normalsGeneric->GetTuple(surroundingArray[i], normalOfPt);

    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(surroundingArray[i], &getPt[0]);

    // TODO : Need to figure out thickness formula considering pack and
    // surrounding mesh cell size and unit normal scale.
    double newCoords[3];
    newCoords[0] = getPt[0] + 0.001 * normalOfPt[0];
    newCoords[1] = getPt[1] + 0.001 * normalOfPt[1];
    newCoords[2] = getPt[2] + 0.001 * normalOfPt[2];

    dataSetSurr->GetPoints()->SetPoint(surroundingArray[i], newCoords);
  }

  // Pack cells
  for (int i = 0; i < (int)packArray.size(); i++) {
    double normalOfPt[3];
    normalsGeneric->GetTuple(packArray[i], normalOfPt);

    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(packArray[i], &getPt[0]);

    double newCoords[3];
    newCoords[0] = (getPt[0] + 0.001 * normalOfPt[0]);
    newCoords[1] = (getPt[1] + 0.001 * normalOfPt[1]);
    newCoords[2] = (getPt[2] + 0.001 * normalOfPt[2]);

    dataSetSurr->GetPoints()->SetPoint(packArray[i], newCoords);
  }

  // Making cells now and plotting them
  int newCells = (int)(packArray.size() / 4);
  int startNum = (int)(dataSetSurr->GetNumberOfCells());

  // Finally, creating VTK cells
  int it7 = 0;
  std::vector<int> pntCohesiveIds;
  int nCelCohesivePnts = 0;
  for (int icl = 0; icl < newCells; icl++) {
    nCelCohesivePnts = 4;
    pntCohesiveIds.resize((nCelCohesivePnts * 2), -1);
    for (int ip = 0; ip < nCelCohesivePnts; ip++) {
      pntCohesiveIds[ip] = packArray[it7];
      pntCohesiveIds[ip + 4] = surroundingArray[it7];
      it7++;
    }
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(pntCohesiveIds.size());
    for (auto pit = pntCohesiveIds.begin(); pit != pntCohesiveIds.end(); pit++)
      vtkCellIds->SetId(pit - pntCohesiveIds.begin(), *pit);
    dataSetSurr->InsertNextCell(12, vtkCellIds);
  }

  int endNum = dataSetSurr->GetNumberOfCells();

  // ****************************** //
  // Zero volume cells visualization method for debugging purpose.
  std::vector<int> cohCellIds;
  for (int i = startNum; i < endNum; i++) { cohCellIds.push_back(i); }

  int size3 = cohCellIds.size();
  // Takes cell Ids as input and writes VTK mesh
  vtkSmartPointer<vtkUnstructuredGrid> visDataSet =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  std::vector<int> ptIdzz2;
  for (int i = 0; i < size3; i++) {
    vtkCell *cell;
    cell = dataSetSurr->GetCell(cohCellIds[i]);

    vtkIdList *pts = cell->GetPointIds();

    for (int j = 0; j < 8; j++) { ptIdzz2.push_back(pts->GetId(j)); }
  }

  int lvar = 0;

  vtkSmartPointer<vtkPoints> pointsViz = vtkSmartPointer<vtkPoints>::New();

  for (int i = 0; i < (size3 * 8); i++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(ptIdzz2[i], &getPt[0]);
    pointsViz->InsertNextPoint(getPt[0], getPt[1], getPt[2]);
  }

  visDataSet->SetPoints(pointsViz);

  for (int i = 0; i < size3; i++) {
    std::vector<int> ptnewIdz;
    for (int j = 0; j < 8; j++) {
      ptnewIdz.push_back(lvar);
      lvar++;
    }
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(ptnewIdz.size());
    for (auto pit = ptnewIdz.begin(); pit != ptnewIdz.end(); pit++)
      vtkCellIds->SetId(pit - ptnewIdz.begin(), *pit);
    visDataSet->InsertNextCell(12, vtkCellIds);
  }

  auto vm = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
      new NEM::MSH::vtkGeoMesh(visDataSet));
  vm->write("ArtificialElements.vtu");

  // ********************************************* //
  auto vm2 = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
      new NEM::MSH::vtkGeoMesh(dataSetSurr));
  vm2->write(outName);

  // Convert final mesh to tetrahedra
  bool tetra = false;
  std::string ofnameTet = "Tetra" + outName;

  if (tetra) {
    vtkSmartPointer<vtkDataSetTriangleFilter> triFilter =
        vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triFilter->SetInputData(dataSetSurr);
    triFilter->Update();

    auto newData_tet = triFilter->GetOutput();

    auto vm3 = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
        new NEM::MSH::vtkGeoMesh(newData_tet));
    vm3->write(ofnameTet);
  }
}

// Generates Periodic Mesh
void MeshManipulationFoam::periodicMeshMapper(std::string &patch1,
                                              std::string &patch2) {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);

  Foam::fileName one = ".";
  Foam::fileName two = ".";
  Foam::Time runTime(controlDict_.get(), one, two);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // Foam::word regionName = Foam::polyMesh::defaultRegion;

  // auto _fmeshPacks = new Foam::polyMesh(Foam::IOobject(
  //     regionName, runTime.timeName(), runTime, Foam::IOobject::MUST_READ));

  auto fm = std::unique_ptr<NEM::MSH::foamGeoMesh>(
      NEM::MSH::foamGeoMesh::Read(".foam"));

  const auto &_fmeshPacks = fm->getFoamMesh();

  // Creating patch objects for periodic boundary patches
  int patchIDIn = _fmeshPacks.boundaryMesh().findPatchID(patch1);
  int patchIDOut = _fmeshPacks.boundaryMesh().findPatchID(patch2);

  const Foam::polyPatch &InPolyPatch = _fmeshPacks.boundaryMesh()[patchIDIn];
  const Foam::polyPatch &OutPolyPatch = _fmeshPacks.boundaryMesh()[patchIDOut];

  // Point IDs from periodic patches.
  Foam::labelList inPts(InPolyPatch.meshPoints());
  Foam::labelList outPts(OutPolyPatch.meshPoints());

  if (inPts.size() != outPts.size()) {
    std::cerr << "Selected boundaries are not periodic!"
              << "\nExiting!" << std::endl;
    throw;
  }

  Foam::pointField packsPF = _fmeshPacks.points();

  // Creating map of IDs on both patches
  std::map<int, int> periodicMap;

  for (int i = 0; i < inPts.size(); i++) { periodicMap[inPts[i]] = outPts[i]; }

  ofstream myfile;
  myfile.open("PeriodicMap.csv");
  for (int i = 0; i < inPts.size(); i++) {
    myfile << inPts[i] << "," << outPts[i] << std::endl;
  }
  myfile.close();

  // Method test using VTK unstructured writer
  // Very useful for extracting face node orders. Might use this in
  // cohesive elements. Also, Foam::facePointPatch::pointNormals()
  // can be useful for artificial thickness elements
  vtkSmartPointer<vtkUnstructuredGrid> dataSetTest =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> pointsMesh = vtkSmartPointer<vtkPoints>::New();

  for (int i = 0; i < packsPF.size(); i++) {
    pointsMesh->InsertNextPoint(packsPF[i].x(), packsPF[i].y(), packsPF[i].z());
  }

  dataSetTest->SetPoints(pointsMesh);

  forAll (_fmeshPacks.boundaryMesh()[patchIDIn], facei) {
    const Foam::label &faceID =
        _fmeshPacks.boundaryMesh()[patchIDIn].start() + facei;
    std::vector<int> pntIds;
    pntIds.resize(4, -1);
    int g = 0;
    forAll (_fmeshPacks.faces()[faceID], nodei) {
      const Foam::label &nodeID = _fmeshPacks.faces()[faceID][nodei];
      pntIds[g] = nodeID;
      g++;
    }
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(pntIds.size());
    for (auto pit = pntIds.begin(); pit != pntIds.end(); pit++)
      vtkCellIds->SetId(pit - pntIds.begin(), *pit);
    dataSetTest->InsertNextCell(9, vtkCellIds);
  }

  forAll (_fmeshPacks.boundaryMesh()[patchIDOut], facei) {
    const Foam::label &faceID =
        _fmeshPacks.boundaryMesh()[patchIDOut].start() + facei;
    std::vector<int> pntIds;
    pntIds.resize(4, -1);
    int g = 0;
    forAll (_fmeshPacks.faces()[faceID], nodei) {
      const Foam::label &nodeID = _fmeshPacks.faces()[faceID][nodei];
      pntIds[g] = nodeID;
      g++;
    }
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(pntIds.size());
    for (auto pit = pntIds.begin(); pit != pntIds.end(); pit++)
      vtkCellIds->SetId(pit - pntIds.begin(), *pit);
    dataSetTest->InsertNextCell(9, vtkCellIds);
  }

  auto vm = std::unique_ptr<NEM::MSH::vtkGeoMesh>(
      new NEM::MSH::vtkGeoMesh(dataSetTest));
  vm->write("PeriodicMesh.vtu");
}