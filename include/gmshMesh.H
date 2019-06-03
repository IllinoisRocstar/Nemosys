#ifndef NEMOSYS_GMSHMESH_H
#define NEMOSYS_GMSHMESH_H

#include <gmsh/GModel.h>

#include "meshBase.H"

class gmshMesh : public meshBase {
  // --- constructors and destructors
  public:
    // default constructor
    gmshMesh();

    // construct from dataset
    gmshMesh(meshBase* mb);

    // initialize from a file name
    explicit gmshMesh(const std::string &fname);

    ~gmshMesh() override;

  // --- diagnostics
  public:
    void report() const override;

  // --- access
  public:
    // get point with id
    std::vector<double> getPoint(int id) const override;
    // get 3 vecs with x,y and z coords
    std::vector<std::vector<double>> getVertCrds() const override;
    // get cell with id : returns point indices and respective coordinates
    std::map<int, std::vector<double>> getCell(int id) const override;
    // get vector of coords of cell with id
    std::vector<std::vector<double>> getCellVec(int id) const override;
    // get edge lengths of dataSet
    void inspectEdges(const std::string& ofname) const override;
    // extract the surface mesh
    vtkSmartPointer<vtkDataSet> extractSurface() override;
    // get diameter of circumsphere of each cell
    std::vector<double> getCellLengths() const override;
    // get center of a cell
    std::vector<double> getCellCenter(int cellID) const override;
    // get cell type as an integer
    // assumes all elements are the same type
    int getCellType() const override;
    // get connectivities. this is only safe to use if mesh has cells of the same type
    // or you have information on the number of cells of each type and the order in
    // which they appear (for look up in resulting vector)
    std::vector<int> getConnectivities() const override;

  // --- read, write and conversion
  public:
    // intended to be used only for reading gmsh mesh files
    void read(const std::string &fname);    
    // write the mesh to file named fname
    void write(const std::string &fname) const override;
    void
    write(const std::string &fname, double mshFileVersion, bool binary) const;

  //   --- Gmsh data types
  private:
//    GModel _gmshGModel;
};

#endif //NEMOSYS_GMSHMESH_H
