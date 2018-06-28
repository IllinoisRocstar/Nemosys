#ifndef meshSrch_H
#define meshSrch_H

// Nemosys
#include "meshBase.H"

class meshSrch : public meshBase
{

  // constructors and destructors
  public:
    meshSrch(meshBase* mb) : vcl(NULL), meshBase(*mb) 
    {}

    ~meshSrch();

  // cell search methods
  public:
    void FindCellsWithinBounds(std::vector<double> bb, std::vector<double>& ids, bool fulImrsd) const;

  // internal management
  private:
    meshSrch();     
    void buildCellLocator();

  private:
    vtkSmartPointer<vtkCellLocator> vcl;

};
#endif
