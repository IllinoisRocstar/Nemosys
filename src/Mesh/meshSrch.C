// Nemosys
#include "meshSrch.H"

// VTK
#include <vtkCellLocator.h>


//meshSrch::meshSrch(meshBase* mb) : vcl(NULL), meshBase(*mb) 
//{}

void meshSrch::buildCellLocator()
{
  if (!vcl)
  {
    // Create the tree
    vcl=vtkSmartPointer<vtkCellLocator>::New();
    vcl->SetDataSet(dataSet);
    vcl->BuildLocator();
  }
}

void meshSrch::FindCellsWithinBounds(std::vector<double> bb, std::vector<double>& ids, bool fulImrsd) const
{
}

