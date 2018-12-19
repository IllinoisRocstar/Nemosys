// Nemosys
#include "meshSrch.H"

// VTK
#include <vtkCellLocator.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>

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
        std::cout << "Remove " << nr << " cells from the list.\n";
    }
    else
        ids = aids;
}

void meshSrch::FindPntsOnTriSrf(std::vector<double>& crds, std::vector<int>& conn, std::set<int>& ids, double tol)
{
    // create polydata
    vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
    for (int iPnt=0; iPnt< (crds.size() / 3); iPnt++)
      pnts->InsertNextPoint(crds[iPnt*3], crds[iPnt*3+1], crds[iPnt*3+2]);
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    for (int iCel=0; iCel< (conn.size() / 3); iCel++)
    {
        polys->InsertNextCell(3);
        polys->InsertCellPoint(conn[iCel*3]);
        polys->InsertCellPoint(conn[iCel*3+1]);
        polys->InsertCellPoint(conn[iCel*3+2]);
    }
    vtkSmartPointer<vtkPolyData> polyData = 
        vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(pnts);
    polyData->SetPolys(polys);

    // Write the file
    //vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    // vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    //writer->SetFileName("test.vtp");
    //writer->SetInputData(polyData);
    //// Optional - set the mode. The default is binary.
    ////writer->SetDataModeToBinary();
    //writer->SetDataModeToAscii();
    //writer->Write();

    // find nodes residing on the trisurf
    // create cell locator
    vtkSmartPointer<vtkCellLocator> cellLocator = 
        vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(polyData);
    cellLocator->BuildLocator();
    
    double testPoint[3];
    double closestPoint[3];
    double closestPointDist2; 
    vtkIdType cellId;
    int subId;
    for (int iPt=0; iPt<getNumberOfPoints(); iPt++)
    {
      std::vector<double> pnt = getPoint(iPt);
      cellLocator->FindClosestPoint(&pnt[0], closestPoint, cellId, subId, closestPointDist2);
      if (closestPointDist2 < tol)
        ids.insert(iPt+1);
    }
}

void meshSrch::FindPntsOnEdge(std::vector<double>& crds, std::set<int>& ids, double tol)
{
    // create polydata
    vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
    for (int iPnt=0; iPnt< (crds.size() / 3); iPnt++)
      pnts->InsertNextPoint(crds[iPnt*3], crds[iPnt*3+1], crds[iPnt*3+2]);
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    polys->InsertNextCell(2);
    polys->InsertCellPoint(0);
    polys->InsertCellPoint(1);    
    vtkSmartPointer<vtkPolyData> polyData = 
        vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(pnts);
    polyData->SetPolys(polys);

    //// Write the file
    //vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    // vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    //writer->SetFileName("edge.vtp");
    //writer->SetInputData(polyData);
    //// Optional - set the mode. The default is binary.
    ////writer->SetDataModeToBinary();
    //writer->SetDataModeToAscii();
    //writer->Write();

    // find nodes residing on the edge
    // create cell locator
    vtkSmartPointer<vtkCellLocator> cellLocator = 
        vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(polyData);
    cellLocator->BuildLocator();
    
    double testPoint[3];
    double closestPoint[3];
    double closestPointDist2; 
    vtkIdType cellId;
    int subId;
    for (int iPt=0; iPt<getNumberOfPoints(); iPt++)
    {
      std::vector<double> pnt = getPoint(iPt);
      cellLocator->FindClosestPoint(&pnt[0], closestPoint, cellId, subId, closestPointDist2);
      if (closestPointDist2 < tol)
        ids.insert(iPt+1);
    }

}


// checks for duplicate elements
bool meshSrch::chkDuplElm() const
{
    std::set<std::vector<int>> ids;
    for (int ic=0; ic<getNumberOfCells(); ic++)
    {
        std::vector<int> cid;
        vtkIdList* idl = vtkIdList::New(); 
        idl = dataSet->GetCell(ic)->GetPointIds();
        for (int id=0; id<idl->GetNumberOfIds(); id++)
            cid.push_back(idl->GetId(id));
        std::pair<std::set<std::vector<int>>::iterator, bool> ret;
        ret = ids.insert(cid);
        if (!ret.second)
        {
            return(true);
        }
    }
    return(false);
}

