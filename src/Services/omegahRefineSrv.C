#include "omegahRefineSrv.H"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <Omega_h_adapt.hpp>
#include <Omega_h_class.hpp>

#include "oshGeoMesh.H"

namespace NEM {
namespace SRV {

Omega_h::Verbosity ParseVerbosity(const std::string &verbosity) {
  if (verbosity == "Silent")
    return Omega_h::SILENT;
  else if (verbosity == "Each Adapt")
    return Omega_h::EACH_ADAPT;
  else if (verbosity == "Each Rebuild")
    return Omega_h::EACH_REBUILD;
  else if (verbosity == "Extra Stats")
    return Omega_h::EXTRA_STATS;

  std::cerr << verbosity << " is not a recognized verbosity level" << std::endl;
  exit(1);
}

Omega_h::VarCompareOpts ParseVarCompareOpts(const std::string &type,
                                            Omega_h::Real tolerance,
                                            Omega_h::Real floor) {
  if (type == "None")
    return {Omega_h::VarCompareOpts::NONE, tolerance, floor};
  else if (type == "Relative")
    return {Omega_h::VarCompareOpts::RELATIVE, tolerance, floor};
  else if (type == "Absolute")
    return {Omega_h::VarCompareOpts::ABSOLUTE, tolerance, floor};

  std::cerr << type << " is not a recognized variable compare type"
            << std::endl;
  exit(1);
}

vtkStandardNewMacro(omegahRefineSrv)

omegahRefineSrv::omegahRefineSrv() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  std::cout << "omegahRefineSrv constructed" << std::endl;
}

omegahRefineSrv::~omegahRefineSrv() {
  std::cout << "omegahRefineSrv destructed" << std::endl;
}

int omegahRefineSrv::FillInputPortInformation(int vtkNotUsed(port),
                                              vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "oshGeoMesh");
  return 1;
}

int omegahRefineSrv::FillOutputPortInformation(int vtkNotUsed(port),
                                               vtkInformation *info) {
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "oshGeoMesh");
  return 1;
}

int omegahRefineSrv::RequestData(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output. Input may just have the geoMeshBase
  // interface, but output should be an Omega_h mesh.
  MSH::oshGeoMesh *input =
      MSH::oshGeoMesh::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  MSH::oshGeoMesh *output =
      MSH::oshGeoMesh::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  Omega_h::Mesh oshOutput = input->getOshMesh();

  if (!input->getOshMesh().is_valid()) return 0;

  // Attempt to classify points and cells, in case not called previously
  Omega_h::finalize_classification(&oshOutput);

  // Define the metric input options
  Omega_h::MetricInput metricInput;
  metricInput.verbose = this->Verbose;
  for (const auto &ms : this->MetricSources)
    metricInput.add_source(
        {ms.type, ms.knob, ms.tag_name, ms.isotropy, ms.scales});
  metricInput.should_limit_lengths = this->ShouldLimitLengths;
  if (this->MaxLength > 0.0) metricInput.max_length = this->MaxLength;
  if (this->MinLength > 0.0) metricInput.min_length = this->MinLength;
  metricInput.should_limit_gradation = this->ShouldLimitGradation;
  if (this->MaxGradationRate > 0.0)
    metricInput.max_gradation_rate = this->MaxGradationRate;
  if (this->GradationConvergenceTolerance > 0.0)
    metricInput.gradation_convergence_tolerance =
        this->GradationConvergenceTolerance;
  metricInput.should_limit_element_count = this->ShouldLimitElementCount;
  if (this->MaxElementCount > 0.0)
    metricInput.max_element_count = this->MaxElementCount;
  if (this->MinElementCount > 0.0)
    metricInput.min_element_count = this->MinElementCount;
  if (this->ElementCountOverRelaxation > 0.0)
    metricInput.element_count_over_relaxation =
        this->ElementCountOverRelaxation;
  if (this->NsmoothingSteps > 0)
    metricInput.nsmoothing_steps = this->NsmoothingSteps;

  // Add implied metric and target metric
  if (!this->MetricSources.empty())
    Omega_h::generate_target_metric_tag(&oshOutput, metricInput);
  Omega_h::add_implied_metric_tag(&oshOutput);

  // Define the adaptation options
  Omega_h::AdaptOpts adaptOpts(&oshOutput);
  if (this->MinLengthDesired > 0.0)
    adaptOpts.min_length_desired = this->MinLengthDesired;
  if (this->MaxLengthDesired > 0.0)
    adaptOpts.max_length_desired = this->MaxLengthDesired;
  if (this->MaxLengthAllowed > 0.0)
    adaptOpts.max_length_allowed = this->MaxLengthAllowed;
  if (this->MinQualityAllowed > 0.0)
    adaptOpts.min_quality_allowed = this->MinQualityAllowed;
  if (this->MinQualityDesired > 0.0)
    adaptOpts.min_quality_desired = this->MinQualityDesired;
  if (this->NsliverLayers > 0) adaptOpts.nsliver_layers = this->NsliverLayers;
  adaptOpts.verbosity = ParseVerbosity(this->Verbosity);
  if (this->LengthHistogramMin > 0.0)
    adaptOpts.length_histogram_min = this->LengthHistogramMin;
  if (this->LengthHistogramMax > 0.0)
    adaptOpts.length_histogram_max = this->LengthHistogramMax;
  if (this->NlengthHistogramBins > 0)
    adaptOpts.nlength_histogram_bins = this->NlengthHistogramBins;
  if (this->NqualityHistogramBins > 0)
    adaptOpts.nquality_histogram_bins = this->NqualityHistogramBins;
  adaptOpts.should_refine = this->ShouldRefine;
  adaptOpts.should_coarsen = this->ShouldCoarsen;
  adaptOpts.should_swap = this->ShouldSwap;
  adaptOpts.should_coarsen_slivers = this->ShouldCoarsenSlivers;
  adaptOpts.should_prevent_coarsen_flip = this->ShouldPreventCoarsenFlip;

  if (TransferOptsTypeMap.empty()) {
    std::map<std::string, Omega_h_Transfer> transferMap;
    const std::set<std::string> reservedTagsPoints{
        "class_id", "class_dim", "momentum_velocity_fixed", "coordinates",
        "warp",     "metric",    "target_metric",           "global"};
    for (Omega_h::Int i = 0; i < oshOutput.ntags(0); ++i) {
      auto tag = oshOutput.get_tag(0, i);
      if (reservedTagsPoints.find(tag->name()) == reservedTagsPoints.end()) {
        bool inherit = true;
        for (Omega_h::Int j = 1; j < oshOutput.dim(); ++j) {
          if (!oshOutput.has_tag(j, tag->name())) {
            inherit = false;
          }
        }
        if (inherit) {
          transferMap[tag->name()] = OMEGA_H_INHERIT;
        } else if (tag->type() == OMEGA_H_REAL) {
          transferMap[tag->name()] = OMEGA_H_LINEAR_INTERP;
        }
      }
    }
    const std::set<std::string> reservedTagsCells{"class_id", "class_dim",
                                                  "size", "quality"};
    for (Omega_h::Int i = 0; i < oshOutput.ntags(oshOutput.dim()); ++i) {
      auto tag = oshOutput.get_tag(oshOutput.dim(), i);
      if (tag->type() == OMEGA_H_REAL &&
          reservedTagsCells.find(tag->name()) == reservedTagsCells.end()) {
        auto it = transferMap.find(tag->name());
        if (it == transferMap.end()) {
          transferMap[tag->name()] = OMEGA_H_POINTWISE;
        }
      }
    }
    adaptOpts.xfer_opts.type_map = transferMap;
  } else {
    adaptOpts.xfer_opts.type_map = TransferOptsTypeMap;
    adaptOpts.xfer_opts.integral_map = TransferOptsIntegralMap;
    for (const auto &idm : this->TransferOptsIntegralDiffuseMap)
      adaptOpts.xfer_opts.integral_diffuse_map[idm.first] = ParseVarCompareOpts(
          idm.second.type, idm.second.tolerance, idm.second.floor);
  }

  // Omega_h::vtk::write_vtu("refine_test_" + std::to_string(0) + ".vtu",
  //                         &oshOutput);

  // int i = 0;
  do {
    Omega_h::adapt(&oshOutput, adaptOpts);
    // ++i;
    // Omega_h::vtk::write_vtu("refine_test_" + std::to_string(i) + ".vtu",
    //                         &oshOutput);
  } while (Omega_h::approach_metric(&oshOutput, adaptOpts));

  output->setOshMesh(&oshOutput);

  return 1;
}

}  // namespace SRV
}  // namespace NEM
