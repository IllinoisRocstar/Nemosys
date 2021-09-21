#ifndef NEMOSYS_FOAMMESH_H_
#define NEMOSYS_FOAMMESH_H_

#ifdef HAVE_CFMSH

#include "Mesh/meshBase.H"

// OpenFOAM
#include <argList.H>
#include <fvOptions.H>

namespace FOAM {

/**
    @brief foamMesh class contains methods for reading and writing mesh in
          OpenFOAM polyMesh format. It also overrides certain methods of
          meshBase for decomposing OpenFOAM polyMesh into VTK unstructured
          database to facilitate convertion to and from OpenFOAM mesh to
          different mesh formats through ConversionDriver.
**/
class foamMesh : public meshBase {
  // --- constructor and deconstructor
 public:
  /**
      @brief foamMesh default constructor. Initializes OpenFOAM environment.
      @param readDB If true, reads mesh database from current location.
                    (default is false)
  **/
  explicit foamMesh(bool readDB = false);

  /** @brief foamMesh alternate constructor. Assigns incoming meshBase mesh
            to internal pointer for manipulation. Creates necessary foam
            dictionaries for runtime environment.
      @param fullMesh meshBase shared pointer for mesh database.
  **/
  foamMesh(std::shared_ptr<meshBase> fullMesh);

  /** @brief foamMesh standard destructor
   **/
  ~foamMesh() override;

  // --- diagnostics
 public:
  // generate a report of the mesh
  void report() const override;

  // --- access
 public:
  /**
      @brief general purpose mesh information read method
      @param fname File name to read the mesh from.
  **/
  void read(const std::string &fname) override;

  /** @brief Reads fvMesh for adaptive mesh refinement process
      @param runTime OpenFOAM Time class object
  **/
  void readAMR(const Foam::Time &runTime);

  /**
      @brief get point with id
      @param id VTK Point ID
      @return Vector of point coordinates <X, Y, Z>
  **/
  std::vector<double> getPoint(nemId_t id) const override;

  /**
      @brief get 3 vecs with x,y and z coords
      @returns Vector with all vertices coordinates in mesh
              <vertices><X,Y,Z>
  **/
  std::vector<std::vector<double>> getVertCrds() const override;

  /**
      @brief get cell with id : returns point indices and respective
            coordinates
      @param id VTK cell ID
      @return a map of cell ID to its bounding point IDs.
  **/
  std::map<nemId_t, std::vector<double>> getCell(nemId_t id) const override;

  /**
      @brief get cell points ids and their coordinates
      @param id VTK cell ID
      @return Vector of cell point ids with their coordinates.
  **/
  std::vector<std::vector<double>> getCellVec(nemId_t id) const override;

  /**
      @brief get diameter of circumsphere of each cell
      @return Vector of diameter of circumsphere of each cell
  **/
  std::vector<double> getCellLengths() const override;

  /**
      @brief get center of a cell
      @param cellID VTK cell ID
      @return Vector of cell center coordinates.
  **/
  std::vector<double> getCellCenter(nemId_t cellID) const override;

  /**
      @brief get cell type as an integer. assumes all elements are the
            same type
  **/
  int getCellType() const override;

  /**
      @brief get edge lengths of dataSet
      @param ofname Output file name
  **/
  void inspectEdges(const std::string &ofname) const override;

  /**
      @brief get vtk cell connectivities in terms of point ids
      @return Vector of points
  **/
  std::vector<nemId_t> getConnectivities() const override;

  using meshBase::write;
  /**
      @brief write the mesh to file named fname
      @param fname Output mesh file name (filename.extension)
  **/
  void write(const std::string &fname) const override;

  // --- processing
 public:
  /**
      @brief Extracts foam mesh surface using VTK data set.
  **/
  vtkSmartPointer<vtkDataSet> extractSurface() override;

  /** @brief Return a polyMesh Object
   **/
  // Foam::polyMesh& getMsh();

  // --- internal methods
 private:
  /**
      @brief Decomposes OpenFOAM mesh into VTK unstructured dataset.
  **/
  void genMshDB();

  /**
      @brief Created necessary dictionaries for OpenFOAM runtime environment.
  **/
  void createFoamDicts();

  // --- openfoam data structure
 private:
  Foam::argList *_args;
  Foam::Time *_runTime;
  Foam::fvMesh *_fmesh;

  // --- internal data
 private:
  /** @brief Shared meshBase pointer for foamMesh constructor and some methods
   **/
  std::shared_ptr<meshBase> _volMB;
};

}  // namespace FOAM

#endif

#endif  // NEMOSYS_FOAMMESH_H_
