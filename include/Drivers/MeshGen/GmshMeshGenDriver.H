#ifndef NEMOSYS_GMSHMESHGENDRIVER_H_
#define NEMOSYS_GMSHMESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/MeshGen/MeshGenDriver.H"
#include "gmshParams.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT GmshMeshGenDriver : public MeshGenDriver {
 public:
  using Files = MeshGenFiles;

  GmshMeshGenDriver(Files files, NEM::GEN::gmshParams params);

  const Files &getFiles() const;
  void setFiles(Files files);
  const NEM::GEN::gmshParams &getParams() const;
  void setParams(NEM::GEN::gmshParams params);
  void execute() const override;

  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(NEM::GEN::gmshParams params);
    NEM::GEN::gmshParams params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshEngine = "gmsh";
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  GmshMeshGenDriver();

  Files files_;
  Opts opts_;

  const Opts &getOpts() const;
  void setOpts(Opts opts);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_GMSHMESHGENDRIVER_H_
