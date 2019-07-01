#include <foamMesh.H>
#include <AuxiliaryFunctions.H>

// openfoam headers
#include <fvCFD.H>
#include <fvMesh.H>
#include <vtkTopo.H>
#include <fileName.H>

// vtk
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellTypes.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCell.h>
#include <vtkExtractEdges.h>
#include <vtkGenericCell.h>
#include <vtkCellIterator.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkTriangleFilter.h>

using namespace FOAM;
using namespace nemAux;

//////////////////////////////////
// foamMesh class 
//////////////////////////////////

foamMesh::foamMesh(bool readDB)
{
    // initialize openfoam
    int argc = 1;
    char** argv = new char*[2];
    argv[0] = new char[100];
    strcpy(argv[0], "NONE");
    _args = new Foam::argList(argc, argv);
    Foam::Info<< "Create time\n" << Foam::endl;
    _runTime = new Foam::Time(Foam::Time::controlDictName, *_args);
    Foam::argList::noParallel();

    // reading mesh database from current location
    if (readDB)
        read("");
}

foamMesh::~foamMesh()
{
    if (_args)
    {
        delete _args;
        _args = NULL;
    }
    if(_runTime)
    {
        delete _runTime;
        _runTime = NULL;
    }
    if(_fmesh)
    {
        delete _fmesh;
        _fmesh = NULL;
    }
    std::cout << "foamMesh destroyed" << std::endl;
}

void foamMesh::read(const std::string& fname)
{
    // openFOAM's mesh information is distributed through
    // multiple files therefore fname is meaningless
    // this implementation reads current working directory
    // and finds mesh data files automatically if they exist
    // TODO: Error out when information are not existing

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
      
      if (fname == "PackMesh")
      {
        regionName = "domain1";
      }
      else if (fname == "SurroundingMesh")
      {
        regionName = "domain0";
      }
      else
      {
        regionName = Foam::fvMesh::defaultRegion;
        Foam::Info
            << "Create mesh for time = "
            << _runTime->timeName() << Foam::nl << Foam::endl;
      }

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

void foamMesh::genMshDB()
{

    // declare vtk dataset
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp 
        = vtkSmartPointer<vtkUnstructuredGrid>::New();


    // creating equivalent vtk topology from fvMesh
    // by default polyhedral cells will be decomposed to 
    // tets and pyramids. Additional points will be added
    // to underlying fvMesh.
    std::cout << "Performing topological decomposition.\n";
    Foam::vtkTopo topo(*_fmesh);

    // point data
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

    // cell data
    std::vector<int> pntIds;
    int nCelPnts = 0;
    for (int icl=0; icl<topo.vertLabels().size(); icl++)
    {
        nCelPnts = topo.vertLabels()[icl].size();
        pntIds.resize(nCelPnts, -1);
        for (int ip=0; ip< nCelPnts; ip++)
            pntIds[ip] = topo.vertLabels()[icl][ip];
        createVtkCell(dataSet_tmp, topo.cellTypes()[icl], pntIds);
    }
    dataSet = dataSet_tmp;

    // Obtaining the mesh data and generating requested output file
    numPoints = dataSet->GetNumberOfPoints();
    numCells = dataSet->GetNumberOfCells();
    hasSizeField = false;
    order = 1;
    std::cout << "Number of points " << numPoints << std::endl;
    std::cout << "Number of cells "<< numCells << std::endl;
}


void foamMesh::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
        const int cellType, std::vector<int>& vrtIds)
{
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(vrtIds.size());
    for (auto pit = vrtIds.begin(); pit!= vrtIds.end(); pit++)
        vtkCellIds->SetId(pit-vrtIds.begin(), *pit);
    dataSet->InsertNextCell(cellType,vtkCellIds);
}


// get point with id
std::vector<double> foamMesh::getPoint(nemId_t id) const
{
  double coords[3];
  dataSet->GetPoint(id, coords);
  std::vector<double> result(coords, coords+3);
  return result;
}

// get 3 vectors with x,y and z coords
std::vector<std::vector<double>> foamMesh::getVertCrds() const
{
  std::vector<std::vector<double>> comp_crds(3);
  for (int i = 0; i < 3; ++i)
  {
    comp_crds[i].resize(numPoints);
  }
  double coords[3];
  for (int i = 0; i < numPoints; ++i)
  {
    dataSet->GetPoint(i,coords);
    comp_crds[0][i] = coords[0];
    comp_crds[1][i] = coords[1];
    comp_crds[2][i] = coords[2];
  }
  return comp_crds;
}

// get cell with id : returns point indices and respective coordinates
std::map<nemId_t, std::vector<double>> foamMesh::getCell(nemId_t id) const
{
  if (id < numCells) 
  {
    std::map<nemId_t,std::vector<double>> cell;
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    point_ids = dataSet->GetCell(id)->GetPointIds();
    vtkIdType num_ids = point_ids->GetNumberOfIds();
    for (vtkIdType i = 0; i < num_ids; ++i)
    {
      nemId_t pntId = point_ids->GetId(i);
      std::vector<double> coord = getPoint(pntId);
      cell.insert(std::pair<nemId_t,std::vector<double>> (pntId,coord));
    }

    return cell;
  }
  else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(1);
  }
}

std::vector<std::vector<double>> foamMesh::getCellVec(nemId_t id) const
{
  if (id < numCells) 
  {
    std::vector<std::vector<double>> cell;
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    point_ids = dataSet->GetCell(id)->GetPointIds();
    vtkIdType num_ids = point_ids->GetNumberOfIds();
    cell.resize(num_ids);
    for (vtkIdType i = 0; i < num_ids; ++i)
    {
      nemId_t pntId = point_ids->GetId(i);
      cell[i] = getPoint(pntId);
    }
    return cell;
  }
  else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(1);
  }

}

void foamMesh::inspectEdges(const std::string& ofname) const
{
  std::ofstream outputStream(ofname);
  if (!outputStream.good())
  {
    std::cerr << "error opening " << ofname << std::endl;
    exit(1);
  }

  vtkSmartPointer<vtkExtractEdges> extractEdges 
    = vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges->SetInputData(dataSet);
  extractEdges->Update();
  
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < extractEdges->GetOutput()->GetNumberOfCells(); ++i)
  {
    extractEdges->GetOutput()->GetCell(i, genCell);
    vtkPoints* points = genCell->GetPoints();
    double p1[3], p2[3];
    points->GetPoint(0,p1);
    points->GetPoint(1,p2);
    double len = std::sqrt(pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2));
    outputStream << len << std::endl;
  } 
}

std::vector<nemId_t> foamMesh::getConnectivities() const
{
  std::vector<nemId_t> connectivities;
  vtkSmartPointer<vtkCellIterator> it 
    = vtkSmartPointer<vtkCellIterator>::Take(dataSet->NewCellIterator()); 
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
  {
    vtkSmartPointer<vtkIdList> pointIds = it->GetPointIds();
    for (vtkIdType i = 0; i < pointIds->GetNumberOfIds(); ++i)
    {
      connectivities.push_back(pointIds->GetId(i));
    } 
  }
  return connectivities;
}

void foamMesh::report() const
{
  if (!dataSet)
  {
    std::cout << "dataSet has not been populated!" << std::endl;
    exit(1);
  }

  typedef std::map<int,int> CellContainer;
  // Generate a report
  std::cout << "Processing the dataset generated from " << filename << std::endl
     << " dataSet contains a " 
     << dataSet->GetClassName()
     <<  " that has " << numCells << " cells"
     << " and " << numPoints << " points." << std::endl;

  CellContainer cellMap;
  for (int i = 0; i < numCells; i++)
  {
    cellMap[dataSet->GetCellType(i)]++;
  }

  CellContainer::const_iterator it = cellMap.begin();
  while (it != cellMap.end())
  {
    std::cout << "\tCell type "
              << vtkCellTypes::GetClassNameFromTypeId(it->first)
              << " occurs " << it->second << " times." << std::endl;
    ++it;
  }

  // Now check for point data
  vtkPointData *pd = dataSet->GetPointData();
  if (pd)
  {
    std::cout << " contains point data with "
         << pd->GetNumberOfArrays()
         << " arrays." << std::endl;
    for (int i = 0; i < pd->GetNumberOfArrays(); i++)
    {
      std::cout << "\tArray " << i << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL") ;
      vtkDataArray* da = pd->GetArray(i);
      std::cout << " with " << da->GetNumberOfTuples() 
                << " values. " << std::endl;
    }
  }
  // Now check for cell data
  vtkCellData *cd = dataSet->GetCellData();
  if (cd)
  {
    std::cout << " contains cell data with " << cd->GetNumberOfArrays()
              << " arrays." << std::endl;
    for (int i = 0; i < cd->GetNumberOfArrays(); i++)
    {
      std::cout << "\tArray " << i << " is named "
                << (cd->GetArrayName(i) ? cd->GetArrayName(i) : "NULL") ;
      vtkDataArray* da = cd->GetArray(i);
      std::cout << " with " << da->GetNumberOfTuples() 
                << " values. " << std::endl;
    }
  }
  // Now check for field data
  if (dataSet->GetFieldData())
  {
    std::cout << " contains field data with "
              << dataSet->GetFieldData()->GetNumberOfArrays()
              << " arrays." << std::endl;
    for (int i = 0; i < dataSet->GetFieldData()->GetNumberOfArrays(); i++)
    {
      std::cout << "\tArray " << i
                << " is named " << dataSet->GetFieldData()->GetArray(i)->GetName();
      vtkDataArray* da = dataSet->GetFieldData()->GetArray(i);
      std::cout << " with " << da->GetNumberOfTuples() 
                << " values. " << std::endl;
    }
  }
}

// get diameter of circumsphere of each cell
std::vector<double> foamMesh::getCellLengths() const
{
  std::vector<double> result;
  result.resize(getNumberOfCells());
  for (nemId_t i = 0; i < getNumberOfCells(); ++i)
  {
    result[i] = std::sqrt(dataSet->GetCell(i)->GetLength2());
  } 
  return result;
}

// get center of a cell
std::vector<double> foamMesh::getCellCenter(nemId_t cellID) const
{
  std::vector<double> center(3,0.0);
  std::vector<std::vector<double>> cell = getCellVec(cellID); 
 
  for (int i = 0; i < cell.size(); ++i)
  {
    center = center + cell[i];
  }
  return (1./cell.size())*center;
}

// returns the cell type
int foamMesh::getCellType() const
{
  return dataSet->GetCellType(0);
}

vtkSmartPointer<vtkDataSet> foamMesh::extractSurface()
{
  // extract surface polygons
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfFilt =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfFilt->SetInputData(dataSet);
  surfFilt->Update();

  // triangulate the surface
  vtkSmartPointer<vtkTriangleFilter> triFilt =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triFilt->SetInputData(surfFilt->GetOutput());
  triFilt->Update();

  return triFilt->GetOutput();
}

// write the mesh to file named fname
void foamMesh::write(const std::string &fname) const
{
    std::cerr << "This method is not impletemented yet!\n";
    throw;
}

