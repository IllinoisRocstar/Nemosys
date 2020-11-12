#ifndef ROCRESTARTDRIVER_H
#define ROCRESTARTDRIVER_H

#include <NemDriver.H>

class meshStitcher;
class cgnsAnalyzer;

namespace NEM {
namespace DRV {

typedef std::pair<std::vector<std::shared_ptr<cgnsAnalyzer>>,
                  std::vector<std::shared_ptr<meshBase>>>
    cgVtPair;
class RocRestartDriver : public NemDriver {
 public:
  RocRestartDriver(const std::vector<std::string> &fluidNamesRm,
                   const std::vector<std::string> &ifluidniNamesRm,
                   const std::vector<std::string> &ifluidnbNamesRm,
                   const std::vector<std::string> &ifluidbNamesRm,
                   const std::vector<std::string> &fluidNamesLts,
                   const std::vector<std::string> &ifluidniNamesLts,
                   const std::vector<std::string> &ifluidnbNamesLts,
                   const std::vector<std::string> &ifluidbNamesLts,
                   const std::vector<std::string> &burnNamesRm,
                   const std::vector<std::string> &iBurnNamesRm,
                   const std::vector<std::string> &burnNamesLts,
                   const std::vector<std::string> &iBurnNamesLts);

  ~RocRestartDriver();
  static RocRestartDriver *readJSON(json inputjson);

 private:
  // remeshed cgns file names
  std::vector<std::string> fluidNamesRm;
  std::vector<std::string> burnNamesRm;
  std::vector<std::string> iBurnNamesRm;
  std::vector<std::string> ifluidniNamesRm;
  std::vector<std::string> ifluidnbNamesRm;
  std::vector<std::string> ifluidbNamesRm;
  // last ts cgns file names
  std::vector<std::string> fluidNamesLts;
  std::vector<std::string> burnNamesLts;
  std::vector<std::string> iBurnNamesLts;
  std::vector<std::string> ifluidniNamesLts;
  std::vector<std::string> ifluidnbNamesLts;
  std::vector<std::string> ifluidbNamesLts;
  // stitches files in order listed above
  std::vector<std::unique_ptr<meshStitcher>> stitchers;
  // stitched last ts meshbase in order listed above
  std::vector<std::shared_ptr<meshBase>> mbObjs;
  // all last ts surface cgns files stitched together (mbObjs[1:-1])
  std::shared_ptr<meshBase> stitchedSurf;
  // remeshed, 0 ts loaded cgns files
  std::vector<std::shared_ptr<cgnsAnalyzer>> fluidRmCg;
  std::vector<std::shared_ptr<cgnsAnalyzer>> burnRmCg;
  std::vector<std::shared_ptr<cgnsAnalyzer>> iBurnRmCg;
  std::vector<std::shared_ptr<cgnsAnalyzer>> ifluidNiRmCg;
  std::vector<std::shared_ptr<cgnsAnalyzer>> ifluidNbRmCg;
  std::vector<std::shared_ptr<cgnsAnalyzer>> ifluidBRmCg;
  // remeshed, 0 ts loaded cgns files converted to vtu
  std::vector<std::shared_ptr<meshBase>> fluidRmMb;
  std::vector<std::shared_ptr<meshBase>> burnRmMb;
  std::vector<std::shared_ptr<meshBase>> iBurnRmMb;
  std::vector<std::shared_ptr<meshBase>> ifluidNiRmMb;
  std::vector<std::shared_ptr<meshBase>> ifluidNbRmMb;
  std::vector<std::shared_ptr<meshBase>> ifluidBRmMb;

  // helpers
 private:
  // populates stitchers
  void stitchCGNS(const std::vector<std::string> &fnames, bool surf);
  // populates cgns and meshbase vectors
  void loadPartCgMb();
  cgVtPair loadCGNS(const std::vector<std::string> &fnames, bool surf);
  // stitch ni, b and nb surfaces
  void stitchSurfaces();
  // transfer stitched data to new meshbase partitions
  void transferStitchedToPartCg(const std::string &transferType);
};

// get cg file names using rocstar convention in glob
inline std::vector<std::string> getCgFNames(const std::string &case_dir,
                                            const std::string &prefix,
                                            const std::string &base_t);

}  // namespace DRV
}  // namespace NEM

#endif
