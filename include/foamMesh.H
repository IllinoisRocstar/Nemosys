#ifndef foam_H
#define foam_H

#ifdef HAVE_CFMSH

#include "meshBase.H"

// vtk
#include <vtkUnstructuredGrid.h>

// OpenFOAM
#include <argList.H>
#include <fvOptions.H>

namespace FOAM {

class foamMesh : public meshBase
{
  public:
    explicit foamMesh(bool readDB = false);
    ~foamMesh() override;

  // --- diagnostics
  public:
    // generate a report of the mesh
    void report() const override;

  // access
  public:
    // general purpose mesh information read method
    void read(const std::string &fname) override;
    // get point with id
    std::vector<double> getPoint(nemId_t id) const override;
    // get 3 vecs with x,y and z coords
    std::vector<std::vector<double>> getVertCrds() const override;
    // get cell with id : returns point indices and respective coordinates
    std::map<nemId_t, std::vector<double>> getCell(nemId_t id) const override;
    std::vector<std::vector<double>> getCellVec(nemId_t id) const override;
    // get diameter of circumsphere of each cell
    std::vector<double> getCellLengths() const override;
    // get center of a cell
    std::vector<double> getCellCenter(nemId_t cellID) const override;
    // get cell type as an integer
    // assumes all elements are the same type
    int getCellType() const override;
    // get edge lengths of dataSet
    void inspectEdges(const std::string &ofname) const override;
    std::vector<nemId_t> getConnectivities() const override;
    // write the mesh to file named fname
    void write(const std::string &fname) const override;

  // processing
  public:
    vtkSmartPointer<vtkDataSet> extractSurface() override;

    // internal methods
  private:
    void genMshDB();
    void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                       int cellType, std::vector<int> &vrtIds);

  // openfoam data structure
  private:
    Foam::argList *_args;
    Foam::Time *_runTime;
    Foam::fvMesh *_fmesh;

  // internal data
  private:
    std::shared_ptr<meshBase> _volMB;
};

} // namespace FOAM

#endif

#endif
