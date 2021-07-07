#ifndef NEMOSYS_NETGENMESHGENDRIVER_H_
#define NEMOSYS_NETGENMESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/MeshGen/MeshGenDriver.H"

#include "MeshGeneration/netgenParams.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT NetgenMeshGenDriver : public MeshGenDriver {
 public:
  using Files = MeshGenFiles;

  NetgenMeshGenDriver(Files files, netgenParams params);

  const Files &getFiles() const;
  void setFiles(Files files);
  const netgenParams &getParams() const;
  void setParams(netgenParams params);
  void execute() const override;

  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(netgenParams params);
    netgenParams params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshEngine = "netgen";
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  NetgenMeshGenDriver();

  Files files_;
  Opts opts_;

  const Opts &getOpts() const;
  void setOpts(Opts opts);
};

}
}

#endif  // NEMOSYS_NETGENMESHGENDRIVER_H_
