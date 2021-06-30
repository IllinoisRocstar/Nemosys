#include "Drivers/Refine/OmegahRefineDriver.H"

#include "geoMeshFactory.H"
#include "omegahRefineSrv.H"
#include "oshGeoMesh.H"

namespace NEM {
namespace DRV {

OmegahRefineDriver::Transfer::Transfer(std::string arrayName,
                                             Omega_h_Transfer method)
    : arrayName(std::move(arrayName)), method(method) {}

OmegahRefineDriver::VarCompare::VarCompare(std::string integralName)
    : integralName(std::move(integralName)) {
  // Default values from Omega_h::VarCompareOpts::defaults()
  this->type = "Relative";
  this->tolerance = 1e-6;
  this->floor = 0.;
}

OmegahRefineDriver::Opts::Opts(
    std::vector<NEM::SRV::omegahRefineMetricSource> sources)
    : MetricSources(std::move(sources)) {}

OmegahRefineDriver::OmegahRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

OmegahRefineDriver::OmegahRefineDriver()
    : OmegahRefineDriver({{}, {}}, Opts{{}}) {}

const OmegahRefineDriver::Opts &OmegahRefineDriver::getOpts() const {
  return opts_;
}

void OmegahRefineDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

void OmegahRefineDriver::execute() const {
  auto inMesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::Read(this->files_.inputMeshFile));
  if (this->opts_.reconstructGeo) {
    inMesh->reconstructGeo();
  }
  vtkSmartPointer<NEM::MSH::oshGeoMesh> osh_mesh =
      vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  osh_mesh->takeGeoMesh(inMesh);

  // Pass through omegahRefineSrv
  vtkSmartPointer<NEM::SRV::omegahRefineSrv> orSrv =
      vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  orSrv->AddInputDataObject(osh_mesh);

  // Set metric input options
  if (this->opts_.Verbose) orSrv->SetVerbose(this->opts_.Verbose.value());
  for (const auto &ms : this->opts_.MetricSources)
    orSrv->AddMetricSource(ms.type, ms.knob, ms.tag_name, ms.isotropy,
                           ms.scales);
  if (this->opts_.ShouldLimitLengths)
    orSrv->SetShouldLimitLengths(this->opts_.ShouldLimitLengths.value());
  if (this->opts_.MaxLength) orSrv->SetMaxLength(this->opts_.MaxLength.value());
  if (this->opts_.MinLength) orSrv->SetMinLength(this->opts_.MinLength.value());
  if (this->opts_.ShouldLimitGradation)
    orSrv->SetShouldLimitGradation(this->opts_.ShouldLimitGradation.value());
  if (this->opts_.MaxGradationRate)
    orSrv->SetMaxGradationRate(this->opts_.MaxGradationRate.value());
  if (this->opts_.GradationConvergenceTolerance)
    orSrv->SetGradationConvergenceTolerance(
        this->opts_.GradationConvergenceTolerance.value());
  if (this->opts_.ShouldLimitElementCount)
    orSrv->SetShouldLimitElementCount(
        this->opts_.ShouldLimitElementCount.value());
  if (this->opts_.MaxElementCount)
    orSrv->SetMaxElementCount(this->opts_.MaxElementCount.value());
  if (this->opts_.MinElementCount)
    orSrv->SetMinElementCount(this->opts_.MinElementCount.value());
  if (this->opts_.ElementCountOverRelaxation)
    orSrv->SetElementCountOverRelaxation(
        this->opts_.ElementCountOverRelaxation.value());
  if (this->opts_.NsmoothingSteps)
    orSrv->SetNsmoothingSteps(this->opts_.NsmoothingSteps.value());

  // Set adapt options
  if (this->opts_.MinLengthDesired)
    orSrv->SetMinLengthDesired(this->opts_.MinLengthDesired.value());
  if (this->opts_.MaxLengthDesired)
    orSrv->SetMaxLengthDesired(this->opts_.MaxLengthDesired.value());
  if (this->opts_.MaxLengthAllowed)
    orSrv->SetMaxLengthAllowed(this->opts_.MaxLengthAllowed.value());
  if (this->opts_.MinQualityAllowed)
    orSrv->SetMinQualityAllowed(this->opts_.MinQualityAllowed.value());
  if (this->opts_.MinQualityDesired)
    orSrv->SetMinQualityDesired(this->opts_.MinQualityDesired.value());
  if (this->opts_.NsliverLayers)
    orSrv->SetNsliverLayers(this->opts_.NsliverLayers.value());
  if (this->opts_.Verbosity) orSrv->SetVerbosity(this->opts_.Verbosity.value());
  if (this->opts_.LengthHistogramMin)
    orSrv->SetLengthHistogramMin(this->opts_.LengthHistogramMin.value());
  if (this->opts_.LengthHistogramMax)
    orSrv->SetLengthHistogramMax(this->opts_.LengthHistogramMax.value());
  if (this->opts_.NlengthHistogramBins)
    orSrv->SetNlengthHistogramBins(this->opts_.NlengthHistogramBins.value());
  if (this->opts_.NqualityHistogramBins)
    orSrv->SetNqualityHistogramBins(this->opts_.NqualityHistogramBins.value());
  if (this->opts_.ShouldRefine)
    orSrv->SetShouldRefine(this->opts_.ShouldRefine.value());
  if (this->opts_.ShouldCoarsen)
    orSrv->SetShouldCoarsen(this->opts_.ShouldCoarsen.value());
  if (this->opts_.ShouldSwap)
    orSrv->SetShouldSwap(this->opts_.ShouldSwap.value());
  if (this->opts_.ShouldCoarsenSlivers)
    orSrv->SetShouldCoarsenSlivers(this->opts_.ShouldCoarsenSlivers.value());
  if (this->opts_.ShouldPreventCoarsenFlip)
    orSrv->SetShouldPreventCoarsenFlip(
        this->opts_.ShouldPreventCoarsenFlip.value());

  // Add transfer options
  for (const auto &xfer : this->opts_.TransferOpts) {
    if (xfer.integralName) {
      orSrv->AddTransferOpts(xfer.arrayName, xfer.method,
                             xfer.integralName.value());
    } else {
      orSrv->AddTransferOpts(xfer.arrayName, xfer.method);
    }
  }
  for (const auto &xfer : this->opts_.TransferOptsIntegralDiffuse)
    orSrv->AddTransferOptsIntegralDiffuse(xfer.integralName, xfer.type,
                                          xfer.tolerance, xfer.floor);

  orSrv->Update();

  // Populate the output mesh
  auto outMesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::New(this->files_.outputMeshFile));
  outMesh->takeGeoMesh(orSrv->GetOutput());
  outMesh->write(this->files_.outputMeshFile);
}

}  // namespace DRV
}  // namespace NEM
