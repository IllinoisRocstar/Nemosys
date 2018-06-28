#ifndef meshSrch_H
#define meshSrch_H

// Nemosys
#include "meshBase.H"

class meshSrch : public meshBase
{

  // constructors and destructors
  public:
    meshSrch(meshBase* mb) : 
    upd_vcl(true), meshBase(*mb)
    {};
    
    static meshSrch* Create(meshBase* mb)     
    {
      meshSrch* ms = new meshSrch(mb); 
      return(ms);
    }

    static std::shared_ptr<meshSrch> CreateShared(meshBase* mb)
    {
      std::shared_ptr<meshSrch> ms;
      ms.reset(meshSrch::Create(mb));
      return ms;
    }

    static std::unique_ptr<meshSrch> CreateUnique(meshBase* mb)
    {
      return std::unique_ptr<meshSrch>(meshSrch::Create(mb));
    }


    ~meshSrch() {};

  // general access
  public:
    std::vector<double> getCellCenter(int cellID) const;

  // cell search methods
  public:
    void FindCellsWithinBounds(std::vector<double>& bb, std::vector<int>& ids, bool fulImrsd=true);

  // misc
  public:
    std::vector<std::vector<double>> getCellVec(int id) const;
    std::vector<double> getPoint(int id) const;    
    int getCellType() const {};

  // internal management
  private:
    meshSrch();     
    void buildCellLocator();

  private:
    bool upd_vcl;
    vtkSmartPointer<vtkCellLocator> vcl;

};
#endif
