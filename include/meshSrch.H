#ifndef meshSrch_H
#define meshSrch_H

// Nemosys headers
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

  // check for special conditions
  public:
    bool chkDuplElm() const; // finds duplicate elements

  // cell search methods
  public:
    void FindCellsWithinBounds(std::vector<double>& bb, std::vector<int>& ids, bool fulImrsd=true);
    // get coordinates and connectivities of the surface triangulation
    // and returns ids for the nodes that reside on the triangulation with
    // given tolerance
    void FindPntsOnTriSrf(std::vector<double>& crds, std::vector<int>& conn, std::set<int>& ids, double tol=0.1e-15);
    // get coordinates of the start and end points of an edge 
    // and returns ids for the nodes that reside on the edge with
    // given tolerance
    void FindPntsOnEdge(std::vector<double>& crds, std::set<int>& ids, double tol=0.1e-15);

  // misc
  public:
    std::vector<std::vector<double>> getCellVec(int id) const override;
    std::vector<double> getPoint(int id) const override;
    int getCellType() const override {};
    std::vector<std::vector<double>> getVertCrds() const override {};
    std::map<int, std::vector<double>> getCell(int id) const override {};
    void inspectEdges(const std::string& ofname) const override {};
    vtkSmartPointer<vtkDataSet> extractSurface() override {};
    int getCellDataIdx(const char *name) override {};
    std::vector<double> getCellLengths() const override {};
    std::vector<int> getConnectivities() const override {};

    void write(const std::string &fname) const override {};

  // internal management
  private:
    meshSrch();     
    void buildCellLocator();

  private:
    bool upd_vcl;
    vtkSmartPointer<vtkCellLocator> vcl;

};
#endif
