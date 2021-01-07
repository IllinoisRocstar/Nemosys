#ifndef GMSHGEN_H
#define GMSHGEN_H

#include "meshGen.H"

namespace NEM {

namespace GEN {

class gmshParams;

class gmshGen : public meshGen {
 public:
  /** @brief gmshGen default constructor
   **/
  gmshGen();

  /** @brief gmshGen alternate constructor with parameters
@param params gmshParams object
**/
  gmshGen(gmshParams *params);

  /** @brief gmshGen standard destructor
   **/
  ~gmshGen();

 public:
  /** @brief Creates mesh from input STEP file
      @param fname The input STEP file
  **/
  int createMeshFromSTL(const char *fname) override;

  int createMeshFromSTEP(const char *fname);

 public:

  /** @brief Get filename of generated .msh file
   **/
  const std::string &getMshFilename() const { return mshFname; }

 private:
  /** @brief Sets the global geometry and meshing options
   **/
  void globalOptions();

  /** @brief Gets geometry entitiy names of STEP geometry
   **/
  void getGeomNames();

  /** @brief Gets the surface colors of STEP geometry
   **/
  void getSurfaceColors();

  /** @brief Applies mesh size fields
   **/
  void meshSizeFields();

  /** @brief Applies physical names based on color
   **/
  void applyColorNames();

  /** @brief Applies transfinite settings to prescribed hexahedral volumes
   **/
  void applyTransfiniteVolumes();

  /** @brief gmshParams object Parameters
   **/
  gmshParams *meshParams;  // Gmsh Parameters

  /** @brief Name of .msh file output by gmsh
   */
  std::string mshFname;
};

}  // namespace GEN

}  // namespace NEM

#endif
