#ifndef NEMOSYS_OMEGAHREFINESRV_H
#define NEMOSYS_OMEGAHREFINESRV_H

#include "srvBase.H"

#include <map>
#include <utility>

#include <Omega_h_defines.hpp>

namespace NEM {
namespace SRV {

/**
 * Store variables to create Omega_h's MetricSource object used to define
 * metrics for refinement
 */
struct omegahRefineMetricSource {
  Omega_h_Source type;   ///< type of source, e.g., by variation, by gradient
  Omega_h::Real knob;    ///< extra parameter depending on source type
  std::string tag_name;  ///< name of field defining the source
  Omega_h_Isotropy isotropy;  ///< controls the conversion of anisotropic
                              ///< metrics to isotropic ones
  Omega_h_Scales scales;      ///< determines whether a metric is allowed to be
                              ///< scaled to satisfy an element count constraint
};

/**
 * Store variables to create Omega_h's VarCompareOpts object used to define
 * integral diffuse map during conservative data transfer
 */
struct omegahRefineVarCompareOpts {
  std::string type;         ///< NONE, RELATIVE, or ABSOLUTE comparison of data
  Omega_h::Real tolerance;  ///< tolerance threshold for comparison
  Omega_h::Real floor;      ///< RELATIVE: values below floor snap to 0.0
};

/**
 * @class omegahRefineSrv
 * @brief A service to use the Omega_h library for refinement
 */
class NEMOSYS_EXPORT omegahRefineSrv : public srvBase {
 public:
  static omegahRefineSrv *New();
  vtkTypeMacro(omegahRefineSrv, srvBase)

  omegahRefineSrv();
  ~omegahRefineSrv() override;

  /** @defgroup metric Metric Input
   *  Metric Input Parameters
   *  @{
   */

  //@{
  /**
   * Turn off the verbosity of Omega_h metric generation.
   */
  vtkSetMacro(Verbose, bool);
  vtkGetMacro(Verbose, bool);
  vtkBooleanMacro(Verbose, bool);
  //@}

  //@{
  /**
   * Add a source to the metric input.
   */
  void AddMetricSource(Omega_h_Source type, Omega_h::Real knob = 1.0,
                       std::string tag_name = "",
                       Omega_h_Isotropy isotropy = OMEGA_H_ANISOTROPIC,
                       Omega_h_Scales scales = OMEGA_H_SCALES) {
    omegahRefineMetricSource orms{type, knob, std::move(tag_name), isotropy,
                                  scales};
    MetricSources.push_back(orms);
  }
  //@}

  //@{
  /**
   * Specify if the metrics should limit the edge lengths. Default is off.
   */
  vtkSetMacro(ShouldLimitLengths, bool);
  vtkGetMacro(ShouldLimitLengths, bool);
  vtkBooleanMacro(ShouldLimitLengths, bool);
  //@}

  //@{
  /**
   * Specify the maximum length of the edges allowed by the metrics.
   */
  vtkSetMacro(MaxLength, Omega_h::Real);
  vtkGetMacro(MaxLength, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the minimum length of the edges allowed by the metrics.
   */
  vtkSetMacro(MinLength, Omega_h::Real);
  vtkGetMacro(MinLength, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify if the metrics should limit the gradation rate. Default is off.
   */
  vtkSetMacro(ShouldLimitGradation, bool);
  vtkGetMacro(ShouldLimitGradation, bool);
  vtkBooleanMacro(ShouldLimitGradation, bool);
  //@}

  //@{
  /**
   * Specify the maximum gradation rate allowed by the metrics.
   */
  vtkSetMacro(MaxGradationRate, Omega_h::Real);
  vtkGetMacro(MaxGradationRate, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the tolerance used when limiting the gradation rate in the metrics.
   */
  vtkSetMacro(GradationConvergenceTolerance, Omega_h::Real);
  vtkGetMacro(GradationConvergenceTolerance, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify if the metrics should limit the element count.
   */
  vtkSetMacro(ShouldLimitElementCount, bool);
  vtkGetMacro(ShouldLimitElementCount, bool);
  vtkBooleanMacro(ShouldLimitElementCount, bool);
  //@}

  //@{
  /**
   * Specify the maximum element count allowed by the metrics.
   */
  vtkSetMacro(MaxElementCount, Omega_h::Real);
  vtkGetMacro(MaxElementCount, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the minimum element count allowed by the metrics.
   */
  vtkSetMacro(MinElementCount, Omega_h::Real);
  vtkGetMacro(MinElementCount, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the over relaxation tolerance when limiting the element count.
   */
  vtkSetMacro(ElementCountOverRelaxation, Omega_h::Real);
  vtkGetMacro(ElementCountOverRelaxation, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the number of smoothing steps when generating the metric.
   */
  vtkSetMacro(NsmoothingSteps, Omega_h::Int);
  vtkGetMacro(NsmoothingSteps, Omega_h::Int);
  //@}

  /** @} */  // end of metric

  /** @defgroup adapt Adaptation Options
   *  @{
   */

  //@{
  /**
   * Specify the minimum desired length.
   */
  vtkSetMacro(MinLengthDesired, Omega_h::Real);
  vtkGetMacro(MinLengthDesired, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the maximum desired length.
   */
  vtkSetMacro(MaxLengthDesired, Omega_h::Real);
  vtkGetMacro(MaxLengthDesired, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the maximum allowed length.
   */
  vtkSetMacro(MaxLengthAllowed, Omega_h::Real);
  vtkGetMacro(MaxLengthAllowed, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the minimum allowed length.
   */
  vtkSetMacro(MinQualityAllowed, Omega_h::Real);
  vtkGetMacro(MinQualityAllowed, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the minimum desired quality.
   */
  vtkSetMacro(MinQualityDesired, Omega_h::Real);
  vtkGetMacro(MinQualityDesired, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the number of sliver layers.
   */
  vtkSetMacro(NsliverLayers, Omega_h::Int);
  vtkGetMacro(NsliverLayers, Omega_h::Int);
  //@}

  //@{
  /**
   * Specify the verbosity of the adaptation method. Options are:
   * - Silent
   * - Each Adapt : after each adaptation step
   * - Each Rebuild : also after each rebuild between steps
   * - Extra Stats : all the above and detailed histogram with stats
   */
  vtkSetMacro(Verbosity, const std::string &);
  vtkGetMacro(Verbosity, const std::string &);
  //@}

  //@{
  /**
   * Specify the minimum depth of the length histogram.
   */
  vtkSetMacro(LengthHistogramMin, Omega_h::Real);
  vtkGetMacro(LengthHistogramMin, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the maximum depth of the length histogram.
   */
  vtkSetMacro(LengthHistogramMax, Omega_h::Real);
  vtkGetMacro(LengthHistogramMax, Omega_h::Real);
  //@}

  //@{
  /**
   * Specify the number of bins in the length histogram.
   */
  vtkSetMacro(NlengthHistogramBins, Omega_h::Int);
  vtkGetMacro(NlengthHistogramBins, Omega_h::Int);
  //@}

  //@{
  /**
   * Specify the number of bins in the quality histogram.
   */
  vtkSetMacro(NqualityHistogramBins, Omega_h::Int);
  vtkGetMacro(NqualityHistogramBins, Omega_h::Int);
  //@}

  //@{
  /**
   * Specify if adapt should refine elements. Default is true.
   */
  vtkSetMacro(ShouldRefine, bool);
  vtkGetMacro(ShouldRefine, bool);
  vtkBooleanMacro(ShouldRefine, bool);
  //@}

  //@{
  /**
   * Specify if adapt should coarsen elements. Default is true.
   */
  vtkSetMacro(ShouldCoarsen, bool);
  vtkGetMacro(ShouldCoarsen, bool);
  vtkBooleanMacro(ShouldCoarsen, bool);
  //@}

  //@{
  /**
   * Specify if adapt should swap elements. Default is true.
   */
  vtkSetMacro(ShouldSwap, bool);
  vtkGetMacro(ShouldSwap, bool);
  vtkBooleanMacro(ShouldSwap, bool);
  //@}

  //@{
  /**
   * Specify if adapt should coarsen slivers. Default is true.
   */
  vtkSetMacro(ShouldCoarsenSlivers, bool);
  vtkGetMacro(ShouldCoarsenSlivers, bool);
  vtkBooleanMacro(ShouldCoarsenSlivers, bool);
  //@}

  //@{
  /**
   * Specify if adapt should prevent coarsen flips. Default is false.
   */
  vtkSetMacro(ShouldPreventCoarsenFlip, bool);
  vtkGetMacro(ShouldPreventCoarsenFlip, bool);
  vtkBooleanMacro(ShouldPreventCoarsenFlip, bool);
  //@}

  /** @} */  // end of adapt

  /** @defgroup transfer Transfer Options
   *  @{
   */

  /**
   * Add a field to be transferred.
   *
   * When using "Conserve" transfer, specify the name of the integral.
   * E.g., "density" -> "mass". See Daniel Ibanez's thesis for details on the
   * methods. Note if no transfer opts are added, tags are transferred using
   * OMEGA_H_INHERIT, OMEGA_H_LINEAR_INTERP, or OMEGA_H_POINTWISE depending on
   * dimensions and type.
   */
  void AddTransferOpts(const std::string &name, Omega_h_Transfer method,
                       const std::string &integral_name = "") {
    TransferOptsTypeMap[name] = method;
    if (!integral_name.empty()) TransferOptsIntegralMap[name] = integral_name;
  }
  /**
   * When specifying integral quantities, set tolerances.
   */
  void AddTransferOptsIntegralDiffuse(const std::string &integral_name,
                                      const std::string &type,
                                      Omega_h::Real tolerance,
                                      Omega_h::Real floor) {
    TransferOptsIntegralDiffuseMap[integral_name] = {type, tolerance, floor};
  }

  /** @} */  // end of transfer

 protected:
  int RequestDataObject(vtkInformation *request,
                        vtkInformationVector **inputVector,
                        vtkInformationVector *outputVector) override;
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  // see algorithm for more info
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

 private:
  // Metric Input options
  bool Verbose{true};
  std::vector<omegahRefineMetricSource> MetricSources;
  bool ShouldLimitLengths{false};
  Omega_h::Real MaxLength{-1.0};
  Omega_h::Real MinLength{-1.0};
  bool ShouldLimitGradation{false};
  Omega_h::Real MaxGradationRate{-1.0};
  Omega_h::Real GradationConvergenceTolerance{-1.0};
  bool ShouldLimitElementCount{false};
  Omega_h::Real MaxElementCount{-1.0};
  Omega_h::Real MinElementCount{-1.0};
  Omega_h::Real ElementCountOverRelaxation{-1.0};
  Omega_h::Int NsmoothingSteps{-1};

  // Adapt Opts
  Omega_h::Real MinLengthDesired{-1.0};
  Omega_h::Real MaxLengthDesired{-1.0};
  Omega_h::Real MaxLengthAllowed{-1.0};
  Omega_h::Real MinQualityAllowed{-1.0};
  Omega_h::Real MinQualityDesired{-1.0};
  Omega_h::Int NsliverLayers{-1};
  std::string Verbosity{"Each Rebuild"};
  Omega_h::Real LengthHistogramMin{-1.0};
  Omega_h::Real LengthHistogramMax{-1.0};
  Omega_h::Int NlengthHistogramBins{-1};
  Omega_h::Int NqualityHistogramBins{-1};
  bool ShouldRefine{true};
  bool ShouldCoarsen{true};
  bool ShouldSwap{true};
  bool ShouldCoarsenSlivers{true};
  bool ShouldPreventCoarsenFlip{false};

  // Transfer Opts
  std::map<std::string, Omega_h_Transfer> TransferOptsTypeMap;
  std::map<std::string, std::string> TransferOptsIntegralMap;
  std::map<std::string, omegahRefineVarCompareOpts>
      TransferOptsIntegralDiffuseMap;
};

}  // namespace SRV
}  // namespace NEM

#endif  // NEMOSYS_OMEGAHREFINESRV_H
