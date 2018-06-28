// Nemosys
#include "meshSrch.H"

// VTK
#include <vtkCellLocator.h>

using namespace nemAux;


// get point with id
std::vector<double> meshSrch::getPoint(int id) const
{
  double coords[3];
  dataSet->GetPoint(id, coords);
  std::vector<double> result(coords, coords+3);
  return result;
}

// returns coordinates of the cell vertices in a vector
std::vector<std::vector<double>> meshSrch::getCellVec(int id) const
{
  if (id < numCells) 
  {
    std::vector<std::vector<double>> cell;
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    point_ids = dataSet->GetCell(id)->GetPointIds();
    int num_ids = point_ids->GetNumberOfIds();
    cell.resize(num_ids);
    for (int i = 0; i < num_ids; ++i) 
    {
      int pntId = point_ids->GetId(i);
      cell[i] = getPoint(pntId);
    }
    return cell;
  }
  else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(1);
  }

}



// get center of a cell
std::vector<double> meshSrch::getCellCenter(int cellID) const
{
  std::vector<double> center(3);
  std::vector<std::vector<double>> cell = getCellVec(cellID); 
 
  for (int i = 0; i < cell.size(); ++i)
    center = center + cell[i];
  return (1./cell.size())*center;
}

void meshSrch::buildCellLocator()
{
  if (upd_vcl)
  {
    // Create the tree
    vcl=vtkSmartPointer<vtkCellLocator>::New();
    vcl->SetDataSet(dataSet);
    vcl->BuildLocator();
    upd_vcl = false;
  }
}

void meshSrch::FindCellsWithinBounds(std::vector<double>& bb, std::vector<int>& ids, bool fulImrsd)
{
    // finding all intersecting cells
    buildCellLocator();
    vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
    vcl->FindCellsWithinBounds(&bb[0], idl);
    std::cout << "Found " << idl->GetNumberOfIds() << " cells.\n";
    std::vector<int> aids;
    for (int idx=0; idx<idl->GetNumberOfIds(); idx++)
        aids.push_back(idl->GetId(idx));
    // removing cells centered out of the bounding box
    int nr = 0;
    if (fulImrsd)
    {
        for (auto it=aids.begin(); it!=aids.end(); it++)
        {
            if ( !isInBBox( getCellCenter(*it), bb) )
            {
                nr++;
                continue;
            }
            else
                ids.push_back(*it);
        }
        std::cout << "Removed additional " << nr << std::endl;
    }
    else
        ids = aids;
}

