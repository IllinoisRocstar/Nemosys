#ifndef SURFACESPLITBYTOPOLOGY_H
#define SURFACESPLITBYTOPOLOGY_H

#include "Time.H"
#include "argList.H"
#include "fileName.H"
#include "fvMesh.H"
#include "triSurface.H"

class surfaceSplitByTopology {
 public:
  /** @brief Constructor */
  surfaceSplitByTopology();

  /** @brief Destructor */
  ~surfaceSplitByTopology();

  /**
   * @brief Splits STL into discrete parts based on topology
   * @param surfFileName Input surface file name
   * @param outFileName  Output surface file name
   * @return Number of discrete geometries identified
   */
  int execute(const Foam::fileName surfFileName,
              const Foam::fileName outFileName);
};

#endif  // SURFACESPLITBYTOPOLOGY_H
