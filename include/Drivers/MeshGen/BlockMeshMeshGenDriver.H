#ifndef NEMOSYS_BLOCKMESHMESHGENDRIVER_H_
#define NEMOSYS_BLOCKMESHMESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/MeshGen/MeshGenDriver.H"
#include "MeshGeneration/blockMeshParams.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT BlockMeshMeshGenDriver : public MeshGenDriver {
 public:
  using Files = DriverOutFile;

  BlockMeshMeshGenDriver(Files file, blockMeshParams params);

  const Files &getFiles() const;
  void setFiles(Files file);
  const blockMeshParams &getParams() const;
  void setParams(blockMeshParams params);
  void execute() const override;

  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(blockMeshParams params);
    blockMeshParams params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshEngine = "blockMesh";
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  BlockMeshMeshGenDriver();

  Files file_;
  Opts opts_;

  const Opts &getOpts() const;
  void setOpts(Opts opts);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_BLOCKMESHMESHGENDRIVER_H_
