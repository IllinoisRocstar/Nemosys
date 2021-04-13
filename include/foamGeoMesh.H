#ifndef NEMOSYS_FOAMGEOMESH_H_
#define NEMOSYS_FOAMGEOMESH_H_

#include "geoMeshBase.H"

#include <vtkUnstructuredGrid.h>

#include <argList.H>
#include <fvOptions.H>
#include <fvMesh.H>

namespace NEM {
namespace MSH {

// TODO: --> Read parameters from filename? (directory, field, runtime)
//       --> Read/write volVectorField/volScalarFields

/**
 * @class foamGeoMesh
 * A concrete implementation of @c geoMeshBase representing a mesh in a @c
 * fvMesh
 */
class NEMOSYS_EXPORT foamGeoMesh : public geoMeshBase {
  public:
   static foamGeoMesh *New();
   vtkTypeMacro(foamGeoMesh, geoMeshBase)

  public:
   /**
    * Create a foamMesh from mesh in current directory
    * @param fileName name of file
    * @param phyGrpArrayName name of a volScalarField defining the
    * physical group of each cell
    * @return new foamGeoMesh from mesh
    */
   static foamGeoMesh *Read(const std::string &fileName,
                            const std::string &phyGrpArrayName = std::string());

  public:
   /**
    * Construct a foamGeoMesh with an empty mesh
    */
   foamGeoMesh();
   
   /**
    * Construct a foamGeoMesh from an existing fvMesh
    * @param foamMesh existing foamMesh
    * @param phyGrpArrayName name of a volScalarField defining the physical 
    * group of each cell
    */
   explicit foamGeoMesh(Foam::fvMesh *foamMesh,
                        const std::string &phyGrpArrayName = std::string());

   /**
    * foamGeoMesh deconstructor
    */
   ~foamGeoMesh() override;

  public:
   /**
    * Writes foam mesh into current directory
    * @param fileName Direcotry path where foam mesh should be written
    */
   void write(const std::string &fileName) override;

   /**
    * Reports data about mesh
    */
   void report(std::ostream &out) const override;

  public:
   /**
   * Copy the mesh
   * @return fvMesh
   */
   const Foam::fvMesh &getFoamMesh() const { return *fmesh_; }
  
   /**
    * Set the foamGeoMesh's mesh to @p foamMesh
    * @param foamMesh existing fvMesh;
    */
   void setFoamMesh(std::unique_ptr<Foam::fvMesh> foamMesh);

   /**
    * Resets the class members and initiates an empty foam mesh
    */
   void resetNative() override;

  private:
   /**
    * Initializes arguments and runtime objects
    */
   void InitializeFoam();

   /**
    * Initializes empty fvMesh object
    */
   std::unique_ptr<Foam::fvMesh> initFoamMesh();

  private:
   /**
    * Create a @c GeoMesh from a @c fvMesh
    * @details Assumes gmsh has already been initialized.
    * @param foamMesh fv mesh
    * @param phyGrpArrayName name of cell data array describing physical group of
    * each cell
    * @return @c GeoMesh
    */
   static GeoMesh foam2GM(Foam::fvMesh *foamMesh,
                          const std::string &phyGrpArrayName = std::string());

   /**
    * Create a @c fvMesh from @c geoMeshBase
    * @param geoMesh GeoMesh
    * @param phyGrpArrayName name of cell data array describing physical group of
    * each cell
    * @return @c fvMesh
    */
   static std::unique_ptr<Foam::fvMesh> GM2foam(
        const GeoMesh &geoMesh,
        Foam::Time *runTime_,
        const std::string &phyGrpArrayName = std::string());

  private:
   // ControlDict
   std::unique_ptr<Foam::dictionary> controlDict_;
  
   // fvSchemes
   std::unique_ptr<Foam::dictionary> fvSchemes_;
  
   // fvSolution
   std::unique_ptr<Foam::dictionary> fvSolution_;
   
   // OpenFOAM time object
   std::unique_ptr<Foam::Time> runTime_;

   // OpenFOAM fvMesh Object
   std::unique_ptr<Foam::fvMesh> fmesh_;

   // vtkUnstructuredGrid
   vtkSmartPointer<vtkUnstructuredGrid> dataSet;
};

} // namespace MSH
} // namespace NEM

#endif  // NEMOSYS_FOAMGEOMESH_H_