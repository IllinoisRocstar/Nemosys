#ifndef meshSrch_H
#define meshSrch_H

// Nemosys headers
#include "nemosys_export.h"
#include "meshBase.H"

#include <set>

class NEMOSYS_EXPORT meshSrch : public meshBase
{
  // constructors and destructors
  public:
    explicit meshSrch(meshBase *mb) : upd_vcl(true), meshBase(*mb) {}

    static meshSrch *Create(meshBase *mb)
    {
      auto *ms = new meshSrch(mb);
      return ms;
    }

    static std::shared_ptr<meshSrch> CreateShared(meshBase *mb)
    {
      std::shared_ptr<meshSrch> ms;
      ms.reset(meshSrch::Create(mb));
      return ms;
    }

    static std::unique_ptr<meshSrch> CreateUnique(meshBase *mb)
    {
      return std::unique_ptr<meshSrch>(meshSrch::Create(mb));
    }

    ~meshSrch() override = default;

  // general access
  public:
    std::vector<std::vector<double>> getCellVec(nemId_t id) const override;
    std::vector<double> getPoint(nemId_t id) const override;
    std::vector<double> getCellCenter(nemId_t cellID) const override;

    int getCellType() const override { return 0; }
    std::vector<std::vector<double>> getVertCrds() const override
    { return std::vector<std::vector<double>>(); }
    std::map<nemId_t, std::vector<double>> getCell(nemId_t id) const override
    { return std::map<nemId_t, std::vector<double>>(); }
    int getCellDataIdx(const std::string &name) override { return 0; }
    std::vector<double> getCellLengths() const override
    { return std::vector<double>(); }
    std::vector<nemId_t> getConnectivities() const override
    { return std::vector<nemId_t>(); }

  // check for special conditions
  public:
    bool chkDuplElm() const; // finds duplicate elements

  // cell search methods
  public:
    void FindCellsWithinBounds(std::vector<double> &bb,
                               std::vector<nemId_t> &ids,
                               bool fulImrsd = true);
    // get coordinates and connectivities of the surface triangulation
    // and returns ids for the nodes that reside on the triangulation with
    // given tolerance
    void FindPntsOnTriSrf(const std::vector<double> &crds,
                          const std::vector<nemId_t> &conn,
                          std::set<nemId_t> &ids,
                          double tol = 0.1e-15) const;
    // get coordinates of the start and end points of an edge 
    // and returns ids for the nodes that reside on the edge with
    // given tolerance
    void FindPntsOnEdge(std::vector<double> &crds,
                        std::set<nemId_t> &ids,
                        double tol = 0.1e-15) const;

  // misc
  public:
    void inspectEdges(const std::string &ofname) const override {}
    vtkSmartPointer<vtkDataSet> extractSurface() override
    { return vtkSmartPointer<vtkDataSet>(); }

    void read(const std::string &fname = std::string()) override {}
    void write(const std::string &fname) const override {}

  // internal management
  private:
    meshSrch();
    void buildCellLocator();

  private:
    bool upd_vcl;
    vtkSmartPointer<vtkCellLocator> vcl;
};

#endif
