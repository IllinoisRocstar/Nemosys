#ifndef MERGEMESHES_H
#define MERGEMESHES_H

#include "Time.H"
#include "argList.H"
#include "fileName.H"
#include "fvMesh.H"
#include "mergePolyMesh.H"
#include "processorMeshes.H"
#include "topoSet.H"

class mergeMeshes {
 public:
  /** @brief Constructor
   **/
  mergeMeshes();

  /** @brief Destructor
   **/
  ~mergeMeshes();

  /** @brief Merges two or more meshes into a single database
   **/
  void execute(const std::string& mainCase, const std::string& addCasePath,
               const std::vector<std::string>& addCases, const int& nDomains,
               const int& dirStat);

 private:
  void getRootCase(Foam::fileName& casePath);
};

#endif