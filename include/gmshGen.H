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
  int createMeshFromSTL(const char *fname);

  /** @brief Creates mesh from input STEP file
      @param fname The input STEP file
  **/
  int createMeshFromSTEP(const char *fname);

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

  /** @brief gmshParams object Parameters
   **/
  gmshParams *meshParams;  // Gmsh Parameters
};

}  // namespace GEN

}  // namespace NEM

#endif
