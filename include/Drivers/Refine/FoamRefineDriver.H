#ifndef NEMOSYS_FOAMREFINEDRIVER_H_
#define NEMOSYS_FOAMREFINEDRIVER_H_
#include "nemosys_export.h"
#include "Drivers/Refine/RefineDriver.H"

namespace NEM {
namespace DRV {

struct NEMOSYS_EXPORT FoamRefineOptsBase {
 public:
  int refineInterval{1};
  int maxRefinement{1};
  int nBufferLayers{1};
  int maxCells{500000};
  bool writeFieldData{false};
  bool writeMesh{false};
  bool writeRefHistory{false};
  double timeStep{1.};
  double endTime{1.};
  double startTime{0.};

  JSONCONS_TYPE_TRAITS_FRIEND

 protected:
  // Hide; should use children
  FoamRefineOptsBase() = default;
  static constexpr const char *method = "FV";
};

class NEMOSYS_EXPORT FoamRefineDriver : public RefineDriver {
 public:
  struct NEMOSYS_EXPORT Opts : public FoamRefineOptsBase {
   public:
    enum class Criteria { VALUE, GRADIENT };
    Opts(Criteria refCriteria, std::string inputFieldFile);
    std::string inputFieldFile{};
    double lowerRefineLevel{-1.};
    double upperRefineLevel{-1.};
    double unrefineAbove{-1.};
    double unrefineBelow{-1.};
    Criteria refCriteria{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };

  FoamRefineDriver(Files files, Opts opts);

  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  FoamRefineDriver();

  Opts opts_;
};

#ifdef MLAMR
class NEMOSYS_EXPORT FoamMLRefineDriver : public RefineDriver {
 public:
  struct NEMOSYS_EXPORT Opts : public FoamRefineOptsBase {
   public:
    Opts() = default;
    std::string mlModel{"fdeep_model.json"};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    static constexpr const char *refCriteria = "ml";
  };

  FoamMLRefineDriver(Files files, Opts opts);

  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  FoamMLRefineDriver();

  Opts opts_;
};
#endif

}  // namespace DRV
}  // namespace NEM
#endif  // NEMOSYS_FOAMREFINEDRIVER_H_
