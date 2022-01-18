#ifndef BLOCKMSH_H
#define BLOCKMSH_H

// Foam headers
#include "IOPtrList.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Time.H"
#include "argList.H"
#include "attachPolyTopoChanger.H"
#include "blockMesh.H"
#include "cellSet.H"
#include "cyclicPolyPatch.H"
#include "foamVtkInternalMeshWriter.H"
#include "fvOptions.H"
#include "polyTopoChange.H"
#include "slidingInterface.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
/**
  @brief blockMsh
**/
class blockMsh {
  // --- constructors and destructor
 public:
  /** @brief blockMsh standard constructor
   **/
  blockMsh();  // Default constructor

  /** @brief blockMsh standard desctructor
   **/
  ~blockMsh();  // Class destructor

  // -- meshGen methods override
 public:
  /** @brief Generates mesh
      @param write Whether to write the foam mesh (true) or keep it
             in memory (false)
   **/
  Foam::fvMesh *generate(
      const std::unique_ptr<Foam::Time> &runTime = NULL,
      const std::unique_ptr<Foam::dictionary> &blockMshDict_ = NULL,
      const bool &write = true);
};

#endif
