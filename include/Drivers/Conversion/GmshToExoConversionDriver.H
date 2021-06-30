#ifndef NEMOSYS_GMSHTOEXOCONVERSIONDRIVER_H_
#define NEMOSYS_GMSHTOEXOCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT GmshToExoConversionDriver : public ConversionDriver {
 public:
  using Files = DriverOutFile;

  struct NEMOSYS_EXPORT PostProcTask {
   public:
    explicit PostProcTask(std::string taskFile);
    std::string taskFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    PostProcTask() = default;
  };

  struct NEMOSYS_EXPORT MeshData {
   public:
    explicit MeshData(std::string meshFile);
    std::string meshFileName{};
    std::string meshName{"default"};
    bool usePhys{false};
    bool makeFreeSurfSS{false};
    bool splitTopBotSS{false};
    std::vector<std::string> sideSetNames{};
    std::map<std::string, std::string> elmBlkNames{};
    jsoncons::optional<std::string> addGlobalNodeSet{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    MeshData() = default;
  };

  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(int numMeshes, std::vector<MeshData> meshData, bool needsPostProc);
    int numMeshes{};
    std::vector<MeshData> meshData{};
    bool needsPostProc{};
    int numTasks{0};
    std::vector<PostProcTask> tasks{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *methodName = "GMSH->EXO";
  };

  GmshToExoConversionDriver(Files file, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files file);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  GmshToExoConversionDriver();

  Files file_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_GMSHTOEXOCONVERSIONDRIVER_H_
