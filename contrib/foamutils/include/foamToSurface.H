#ifndef FOAMTOSURFACE_H
#define FOAMTOSURFACE_H

#include "IOdictionary.H"
#include "MeshedSurfaces.H"
#include "Time.H"
#include "argList.H"
#include "fileName.H"
#include "fvMesh.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "timeSelector.H"

class foamToSurface {
 public:
  /** @brief Constructor
   **/
  foamToSurface();

  /** @brief Destructor
   **/
  ~foamToSurface();

  /** @brief Converts foam mesh to surface representation (STL)
      @param outputfile Output surface file name (STL format)
  **/
  void execute(const std::string &outputfile);
};

#endif