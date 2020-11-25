#ifndef NEMOSYS_REMESHDRIVER_H_
#define NEMOSYS_REMESHDRIVER_H_

#include <map>
#include <string>
#include <vector>

#include "nemosys_export.h"
#include "MeshGenDriver.H"
#include "NemDriver.H"
#include "cgnsAnalyzer.H"
#include "meshBase.H"
#include "meshStitcher.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT RemeshDriver : public NemDriver {
 public:
  RemeshDriver(std::vector<std::string> fluidNames,
               std::vector<std::string> ifluidniNames,
               std::vector<std::string> ifluidnbNames,
               std::vector<std::string> ifluidbNames,
               std::vector<std::string> _burnNames,
               std::vector<std::string> _iBurnNames,
               const jsoncons::json &remeshjson, int numPartitions,
               const std::string &base_t, int writeIntermediateFiles,
               double searchTolerance, const std::string &caseName,
               const std::map<std::string, std::vector<int>> &surfacePatchTypes,
               const std::string &_prefix_path = std::string());
  ~RemeshDriver() override;
  static RemeshDriver *readJSON(const jsoncons::json &inputjson);

 private:
  // CGNS file names
  std::vector<std::string> fluidNames;
  std::vector<std::string> ifluidniNames;
  std::vector<std::string> ifluidnbNames;
  std::vector<std::string> ifluidbNames;
  std::vector<std::string> burnNames;
  std::vector<std::string> iBurnNames;
  // mesh generator which uses remeshing parameters from json
  std::unique_ptr<MeshGenDriver> mshgendrvr;
  // stitches files in order listed above
  std::vector<std::unique_ptr<meshStitcher>> stitchers;
  // stitched CGNS in order listed above
  std::vector<std::shared_ptr<cgnsAnalyzer>> cgObjs;
  // stitched meshBase in order listed above
  std::vector<std::shared_ptr<meshBase>> mbObjs;
  // all surface CGNS files stitched together
  std::shared_ptr<meshBase> stitchedSurf;
  // surface burn files stitched together
  std::shared_ptr<meshBase> stitchedBurnSurf;
  // stitched volume CGNS remeshed
  std::shared_ptr<meshBase> remeshedVol;
  // surface extracted from remeshed stitched volume
  std::shared_ptr<meshBase> remeshedSurf;
  // number of partitions
  int numPartitions;
  // prefix path
  std::string prefixPath;

  // helpers
 private:
  // stitch surf/vol CGNS files, using rocstarCgns if surf=1
  void stitchCGNS(const std::vector<std::string> &fnames, bool surf);
  // creates remeshed volume and surface
  void remesh(const jsoncons::json &remeshjson);
  // stitches the stitched CGNS surfaces into one meshBase
  void stitchSurfaces();
};

// get cg file names using rocstar convention in regex glob
inline std::vector<std::string> getCgFNames(const std::string &case_dir,
                                            const std::string &prefix,
                                            const std::string &base_t);

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_REMESHDRIVER_H_
