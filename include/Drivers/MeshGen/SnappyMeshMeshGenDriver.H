#ifndef NEMOSYS_SNAPPYMESHMESHGENDRIVER_H_
#define NEMOSYS_SNAPPYMESHMESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/MeshGen/MeshGenDriver.H"
#include "MeshGeneration/snappymeshParams.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT SnappyMeshMeshGenDriver : public MeshGenDriver {
 public:
  using Files = MeshGenFiles;

  SnappyMeshMeshGenDriver(Files files, snappymeshParams params);

  const Files &getFiles() const;
  void setFiles(Files files);
  const snappymeshParams &getParams() const;
  void setParams(snappymeshParams params);
  void execute() const override;

  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(snappymeshParams params);
    snappymeshParams params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshEngine = "snappyHexMesh";
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  SnappyMeshMeshGenDriver();

  Files files_;
  Opts opts_;

  const Opts &getOpts() const;
  void setOpts(Opts opts);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_SNAPPYMESHMESHGENDRIVER_H_
