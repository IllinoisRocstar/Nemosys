#ifndef NEMOSYS_AMRFOAM_H_
#define NEMOSYS_AMRFOAM_H_

#include <Switch.H>
#include <bitSet.H>
#include <dynamicFvMesh.H>
#include <hexRef8.H>
#include <motionSolver.H>

#include <vector>

namespace Foam {

/** @brief AMRFoam class is a standalone adaptive mesh refinement application
          that uses polyMesh as its base and provides useful interface for
          use of different solvers and mesh formats. Initially the supported
          refinement criterias would be size field defined over mesh as well
          as a gradient field of any simulation quantity.

          Any object/pointer created for this class in C++ runtime will act
          as finiteVolume mesh of OpenFOAM.
**/
// class motionSolver;
class AMRFoam : public dynamicFvMesh {
  // Constructor and destructor
 public:
  /** @brief Standard constructor. Accepts polyMesh object for initialization.
             It needs dynamicMeshDict present in constant folder. Method for
             same is available in MeshManipulationFoam class as
             "createDynamicMeshDict"
      @param iomesh A polyMesh input object
  **/
  AMRFoam(const Foam::IOobject& iomesh);

  /** @brief Class destructor
   **/
  ~AMRFoam() {
    // Nothing
  }

 public:
  /** @brief updateAMR is a master method for this process. It drives the
             refinement and unrefinement procedure between any two given
             timesteps.
      @param refineInterval Refinement interval/frequency
      @param maxRefinement Maximum refinement allowed by user
      @param vFld A scalar field of refinement criteria quantity
      @param lowerRefineLevel Lower refinement level (Based on vFld)
      @param upperRefineLevel Upper refinement level (Based on vFld)
      @param unrefineAbove Criteria for unrefinement above this value
      @param unrefineBelow Criteria for unrefinement below this value
      @param nBufferLayers Buffer layers for slower than 2:1 refinement
      @param maxCells Maximum number of cells allowed in mesh. Once this level
                      is reached, refinement will be stopped.
      @return A boolean indicating mesh change
  **/
  bool updateAMR(const int& refineInterval, const int& maxRefinement,
                 volScalarField& vFld, const double& lowerRefineLevel,
                 const double& upperRefineLevel, const double& unrefineAbove,
                 const double& unrefineBelow, const int& nBufferLayers,
                 const int& maxCells);

  /** @brief New method for machine learning ML. Takes output from ML model
             and performs refinement/unrefinement.
      @param refineInterval Refinement interval/frequency
      @param maxRefinement Maximum refinement allowed by user
      @param nBufferLayers Buffer layers for slower than 2:1 refinement
      @param maxCells Maximum number of cells allowed in mesh. Once this level
                      is reached, refinement will be stopped.
      @param vFld Machine learning model output vector
      @return A boolean indicating mesh change
  **/
  bool updateAMRML(const int& refineInterval, const int& maxRefinement,
                   const int& nBufferLayers, const int& maxCells,
                   volScalarField& vFld);

  /** @brief A mesh motion solver. This method allows solving of mesh motion
             along with mesh refinment. Flux correction takes place after motion
             and refinement operations are performed. Experimental method.
  **/
  void solveMotion();

  /** @brief Converts volume field to point field
      @param vFld Volume Field.
      @return Point Field.
  **/
  Foam::scalarField cell2Pt(const scalarField& vFld);

  /** @brief Writes refinement history in current time directory.
   **/
  void writeData();

  /** @brief Access to current mesh at every time step
      @return A polyMesh object
  **/
  // const polyMesh& getMesh();
  const polyMesh& getMesh();

  /** @brief Ability to write current mesh into current time directory
   **/
  void writeMesh();

  /** @brief Ability to explicitely disable motion
   **/
  void disableMotion();

  /** @brief Converts incoming size-field/solution-field into pointScalarField
             defined over all points of mesh
      @param inName
      @return A pointScalarField with appended incoming values
  **/
  pointScalarField readIncomingPtField(const std::string& inName,
                                       const std::string& outName);

  /** @brief Converts pointScalarField to volScalarField
      @param pntF A point field defined on mesh points
      @param outName Object name for output scalar
      @return A Scalar Field defined on cells
  **/
  volScalarField pF2vF(const pointScalarField& pntF,
                       const std::string& outName);

  /** @brief Reads the refinement field and registers it.
      @param fldName Name of field present in current runTime directory
  **/
  volScalarField readInitialField(const std::string& fldName);

  /** @brief Creates an empty pointScalarField based on current mesh
      @param fldName Name of field to create in current runTime directory
  **/
  pointScalarField initPntSField(const std::string& fldName);

  /** @brief Creates an empty volScalarField based on current mesh
      @param fldName Name of field to create in current runTime directory
  **/
  volScalarField initScalarField(const std::string& fldName);

  /** @brief Creates an empty volScalarField for gradient magnitude
      @param fldName Name of field to create in current runTime directory
  **/
  volScalarField initGradientField(const std::string& fldName);

  /** @brief Converts vector to volScalarField for refinement/unrefinement
      @param vec Incoming vector
      @return volScalarField
  **/
  volScalarField assignToVolScalarField(const std::vector<int>& vec);

  /** @brief Access to updated field at every time step
   **/
  void enableUpdatedField();

  /** @brief Access to refinement history data
   **/
  void enableRefHistoryData();

  /** @brief Access to mesh data
   **/
  void enableMeshWriting();

  /** @brief reads incoming field and converts to volScalarField
      @param inName Input file name
      @param outName Object name for volScalarField
      @return Scalar Field
  **/
  volScalarField readIncomingCellField(const std::string& inName,
                                       const std::string& outName);

  /** @brief Transforms incoming scalar field into refinement levels understood
             by OpenFOAM environement. Appends this information to lower, upper
             and unrefine level vectors for AMR
  **/
  void transformField(volScalarField& inField);

  /** @brief Returns a gradient of scalar field over mesh
             param fldGrd Scalar Field for Gradient Calculation
      @return scalar field of gradient over mesh
  **/
  volScalarField getGradient(const volScalarField& fldGrd);

  /** @brief Vector to store lower refinement levels
   **/
  std::vector<double> lowerRefLvl;

  /** @brief Vector to store upper refinement levels
   **/
  std::vector<double> upperRefLvl;

  /** @brief Vector to store unrefine levels
   **/
  std::vector<double> unrefAboveLvl;

  /** @brief Vector to store unrefine levels
   **/
  std::vector<double> unrefBelowLvl;

  /** @brief Step size for one time refinement
   **/
  int stepSz = 1;

  /**@brief Boolean to check if mesh motion is needed
   **/
  bool enableMotion{false};

  // Private Attributes
 private:
  /** @brief Boolean to enable writing of updated solution field every timestep
   **/
  bool writeField{false};

  /** @brief Boolean to enable writing of refinement history
   **/
  bool writeRefHistory{false};

  /** @brief Boolean to enable writing of mesh
   **/
  bool writeMeshData{false};

  /** @brief Vector to store user defined vector
   **/
  std::vector<double> incomingField;

 protected:
  //- Mesh cutting engine
  hexRef8 meshCutter_;

  //- Dump cellLevel for postprocessing
  bool dumpLevel_;

  //- Fluxes to map
  HashTable<word> correctFluxes_;

  //- Number of refinement/unrefinement steps done so far.
  label nRefinementIterations_;

  //- Protected cells (usually since not hexes)
  bitSet protectedCell_;

  // Protected Member Functions
  /** @brief Checks for availability of mesh motion keywords in dictionary
   **/
  void checkForMotion();

  /** @brief Calculates and stores non-hex cells in protected cells catagory to
             avoid their refinement
  **/
  void calcProtectedCells();

  /** @brief Unrefines selected cells in a mesh.
      @param splitPoints Selected hanging nodes that can be erased to unrefine
             previously refined cells. Selection is performed using criteria and
             protected cells are excluded.
      @return Mapped Polymesh.
  **/
  virtual autoPtr<mapPolyMesh> refine(const labelList&);

  /** @brief Refines selected cells in a mesh.
      @param cellsToRefine Cell list selected for refinement.
      @return Mapped Polymesh.
  **/
  virtual autoPtr<mapPolyMesh> unrefine(const labelList&);

  //- Calculate cells that cannot be refined since would trigger
  //  refinement of protectedCell_ (since 2:1 refinement cascade)
  void calculateProtectedCells(bitSet& unrefineableCell) const;

  //- Read the projection parameters from dictionary
  void readDict();

  //- Calculates approximate value for refinement level so
  //  we don't go above maxCell
  scalar getRefineLevel(const label maxCells, const label maxRefinement,
                        const scalar refineLevel, const scalarField&) const;

  //- Get per cell max of connected point
  scalarField maxPointField(const scalarField&) const;

  //- Get point max of connected cell
  scalarField maxCellField(const volScalarField&) const;

  scalarField cellToPoint(const scalarField& vFld) const;

  scalarField error(const scalarField& fld, const scalar minLevel,
                    const scalar maxLevel) const;

  //- Select candidate cells for refinement
  virtual void selectRefineCandidates(const scalar lowerRefineLevel,
                                      const scalar upperRefineLevel,
                                      const scalarField& vFld,
                                      bitSet& candidateCell) const;

  //- Subset candidate cells for refinement
  virtual labelList selectRefineCells(const label maxCells,
                                      const label maxRefinement,
                                      const bitSet& candidateCell) const;

  /** @brief Marked the cells to unrefine based on user defined levels
      @param unrefineAbove Unrefine cells above this thresold.
      @param unrefineBelow Unrefine cells below this thresold.
      @param markedCell Unrefinable cells list
      @param pFld PointField at current runTime
      @return List of cells/points to unrefine.
  **/
  virtual labelList selectUnrefinePoints(const scalar unrefineAbove,
                                         const scalar unrefineBelow,
                                         const bitSet& markedCell,
                                         const scalarField& pFld) const;

  //- Extend markedCell with cell-face-cell.
  void extendMarkedCells(bitSet& markedCell) const;

  //- Check all cells have 8 anchor points
  void checkEightAnchorPoints(bitSet& protectedCell) const;

  //- Map single non-flux surface<Type>Field
  //  for new internal faces (e.g. AMR refine). This currently
  //  interpolates values from surrounding faces (faces on
  //  neighbouring cells) that do have a value.
  template <class T>
  void mapNewInternalFaces(const labelList& faceMap,
                           GeometricField<T, fvsPatchField, surfaceMesh>&);

  //- Correct surface fields for new faces
  template <class T>
  void mapNewInternalFaces(const labelList& faceMap);

  //- Correct surface fields for new faces. Converts any oriented
  //  fields into non-oriented (i.e. phi -> Uf) before mapping
  template <class T>
  void mapNewInternalFaces(const surfaceVectorField& Sf,
                           const surfaceScalarField& magSf,
                           const labelList& faceMap);

 public:
  virtual bool update();

  // Pointer for motion solver
  autoPtr<motionSolver> motionPtr_;

  // Member Functions

  //- Direct access to the refinement engine
  const hexRef8& meshCutter() const { return meshCutter_; }

  //- Cells which should not be refined/unrefined
  const bitSet& protectedCell() const { return protectedCell_; }

  //- Cells which should not be refined/unrefined
  bitSet& protectedCell() { return protectedCell_; }

  //- Write using given format, version and compression
  virtual bool writeObject(IOstream::streamFormat fmt,
                           IOstream::versionNumber ver,
                           IOstream::compressionType cmp,
                           const bool valid) const;

  //- Map all fields in time using given map.
  virtual void mapFields(const mapPolyMesh& mpm);
};

}  // namespace Foam

#endif  // NEMOSYS_AMRFOAM_H_