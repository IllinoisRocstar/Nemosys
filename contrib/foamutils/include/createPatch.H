#ifndef CREATEPATCH_H
#define CREATEPATCH_H

#include "IOPtrList.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "SortableList.H"
#include "Time.H"
#include "argList.H"
#include "cyclicPolyPatch.H"
#include "faceSet.H"
#include "fileName.H"
#include "fvMesh.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyModifyFace.H"
#include "polyTopoChange.H"
#include "processorMeshes.H"
#include "syncTools.H"
#include "wordReList.H"
#include "wordRes.H"

class createPatch {
 public:
  /** @brief Constructor
   **/
  createPatch();

  /** Constructor **/
  ~createPatch();

  /** @brief Executes the patch creation directed by createPatchDict
      @param dirStat domain number to skip (output from splitMshRegions)
      @param firstPtch Name of the main patch
  **/
  void execute(const int& dirStat, const std::string& firstPtch,
               const std::unique_ptr<Foam::dictionary>& cpDict);

 private:
  void changePatchID(const Foam::polyMesh& mesh, const Foam::label faceID,
                     const Foam::label patchID, Foam::polyTopoChange& meshMod);

  void filterPatches(Foam::polyMesh& mesh,
                     const Foam::HashSet<Foam::word>& addedPatchNames);

  void dumpCyclicMatch(const Foam::fileName& prefix,
                       const Foam::polyMesh& mesh);

  void separateList(const Foam::vectorField& separation,
                    Foam::UList<Foam::vector>& field);

  template <class CombineOp>
  void syncPoints(const Foam::polyMesh& mesh, Foam::pointField& points,
                  const CombineOp& cop, const Foam::point& nullValue);
};

#endif