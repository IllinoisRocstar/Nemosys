#ifndef NEMOSYS_CFMESHMESHGENDRIVER_H_
#define NEMOSYS_CFMESHMESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/MeshGen/MeshGenDriver.H"
#include "cfmeshParams.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT CFMeshMeshGenDriver : public MeshGenDriver {
 public:
  using Files = MeshGenFiles;

  CFMeshMeshGenDriver(Files files, cfmeshParams params);

  const Files &getFiles() const;
  void setFiles(Files files);
  const cfmeshParams &getParams() const;
  void setParams(cfmeshParams params);
  void execute() const override;

  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(cfmeshParams params);
    cfmeshParams params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshEngine = "cfmesh";
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  CFMeshMeshGenDriver();

  Files files_;
  Opts opts_;

  const Opts &getOpts() const;
  void setOpts(Opts opts);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_CFMESHMESHGENDRIVER_H_
