#include "omegahRefineDriver.H"

#include "geoMeshFactory.H"
#include "omegahRefineSrv.H"
#include "oshGeoMesh.H"

namespace NEM {
namespace DRV {

Omega_h_Source Source(const std::string &source) {
  if (source == "Constant")
    return OMEGA_H_CONSTANT;
  else if (source == "Variation")
    return OMEGA_H_VARIATION;
  else if (source == "Derivative")
    return OMEGA_H_DERIVATIVE;
  else if (source == "Given")
    return OMEGA_H_GIVEN;
  else if (source == "Implied")
    return OMEGA_H_IMPLIED;
  else if (source == "Curvature")
    return OMEGA_H_CURVATURE;

  std::cerr << source << " is not a recognized metric source type" << std::endl;
  exit(1);
}

Omega_h_Isotropy Isotropy(const std::string &isotropy) {
  if (isotropy == "Anisotropic")
    return OMEGA_H_ANISOTROPIC;  // keep anisotropy
  else if (isotropy == "Length")
    return OMEGA_H_ISO_LENGTH;  // use smallest length
  else if (isotropy == "Size")
    return OMEGA_H_ISO_SIZE;  // use equivalent volume

  std::cerr << isotropy << " is not a recognized metric source isotropy"
            << std::endl;
  exit(1);
}

/* determines whether a metric is allowed to be scaled
 * to satisfy an element count constraint
 */
Omega_h_Scales Scales(const std::string &scales) {
  if (scales == "Absolute")
    return OMEGA_H_ABSOLUTE;
  else if (scales == "Scales")
    return OMEGA_H_SCALES;

  std::cerr << scales << " is not a recognized metric source scales"
            << std::endl;
  exit(1);
}

Omega_h_Transfer Transfer(const std::string &transfer) {
  if (transfer == "Inherit")
    return OMEGA_H_INHERIT;
  else if (transfer == "Linear Interp")
    return OMEGA_H_LINEAR_INTERP;
  else if (transfer == "Metric")
    return OMEGA_H_METRIC;
  else if (transfer == "Density")
    return OMEGA_H_DENSITY;
  else if (transfer == "Conserve")
    return OMEGA_H_CONSERVE;
  /* The following require the maps in Omega_h::TransferOpts and are handled
   * separately.
  else if (transfer == "Momentum Velocity")
    return OMEGA_H_MOMENTUM_VELOCITY;
  */
  else if (transfer == "Pointwise")
    return OMEGA_H_POINTWISE;

  std::cerr << transfer << " is not a recognized transfer method" << std::endl;
  exit(1);
}

omegahRefineDriver::omegahRefineDriver(NEM::MSH::geoMeshBase *in_mesh,
                                       NEM::MSH::geoMeshBase *out_mesh,
                                       const jsoncons::json &refine_opts) {
  // Convert to oshGeoMesh
  vtkSmartPointer<NEM::MSH::oshGeoMesh> osh_mesh =
      vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  osh_mesh->takeGeoMesh(in_mesh);

  // Pass through omegahRefineSrv
  vtkSmartPointer<NEM::SRV::omegahRefineSrv> orSrv =
      vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  orSrv->AddInputDataObject(osh_mesh);

  // Set metric input options
  if (refine_opts.contains("Verbose"))
    orSrv->SetVerbose(refine_opts["Verbose"].as<bool>());
  if (refine_opts.contains("Metric Sources"))
    for (const auto &ms : refine_opts["Metric Sources"].array_range())
      orSrv->AddMetricSource(
          Source(ms["Type"].as<std::string>()),
          ms.get_with_default<Omega_h::Real>("Knob", 1.0),
          ms.get_with_default<std::string>("Tag Name", ""),
          Isotropy(ms.get_with_default<std::string>("Isotropy", "Anisotropic")),
          Scales(ms.get_with_default<std::string>("Scales", "Scales")));
  if (refine_opts.contains("Should Limit Lengths"))
    orSrv->SetShouldLimitLengths(
        refine_opts["Should Limit Lengths"].as<bool>());
  if (refine_opts.contains("Max Length"))
    orSrv->SetMaxLength(refine_opts["Max Length"].as<Omega_h::Real>());
  if (refine_opts.contains("Min Length"))
    orSrv->SetMinLength(refine_opts["Min Length"].as<Omega_h::Real>());
  if (refine_opts.contains("Should Limit Gradation"))
    orSrv->SetShouldLimitGradation(
        refine_opts["Should Limit Gradation"].as<bool>());
  if (refine_opts.contains("Max Gradation Rate"))
    orSrv->SetMaxGradationRate(
        refine_opts["Max Gradation Rate"].as<Omega_h::Real>());
  if (refine_opts.contains("Gradation Convergence Tolerance"))
    orSrv->SetGradationConvergenceTolerance(
        refine_opts["Gradation Convergence Tolerance"].as<Omega_h::Real>());
  if (refine_opts.contains("Should Limit Element Count"))
    orSrv->SetShouldLimitElementCount(
        refine_opts["Should Limit Element Count"].as<bool>());
  if (refine_opts.contains("Max Element Count"))
    orSrv->SetMaxElementCount(
        refine_opts["Max Element Count"].as<Omega_h::Real>());
  if (refine_opts.contains("Min Element Count"))
    orSrv->SetMinElementCount(
        refine_opts["Min Element Count"].as<Omega_h::Real>());
  if (refine_opts.contains("Element Count Over Relaxation"))
    orSrv->SetElementCountOverRelaxation(
        refine_opts["Element Count Over Relaxation"].as<Omega_h::Real>());
  if (refine_opts.contains("N smoothing Steps"))
    orSrv->SetNsmoothingSteps(
        refine_opts["N smoothing Steps"].as<Omega_h::Int>());

  // Set adapt options
  if (refine_opts.contains("Min Length Desired"))
    orSrv->SetMinLengthDesired(
        refine_opts["Min Length Desired"].as<Omega_h::Real>());
  if (refine_opts.contains("Max Length Desired"))
    orSrv->SetMaxLengthDesired(
        refine_opts["Max Length Desired"].as<Omega_h::Real>());
  if (refine_opts.contains("Max Length Allowed"))
    orSrv->SetMaxLengthAllowed(
        refine_opts["Max Length Allowed"].as<Omega_h::Real>());
  if (refine_opts.contains("Min Quality Allowed"))
    orSrv->SetMinQualityAllowed(
        refine_opts["Min Quality Allowed"].as<Omega_h::Real>());
  if (refine_opts.contains("Min Quality Desired"))
    orSrv->SetMinQualityDesired(
        refine_opts["Min Quality Desired"].as<Omega_h::Real>());
  if (refine_opts.contains("N sliver Layers"))
    orSrv->SetNsliverLayers(refine_opts["N sliver Layers"].as<Omega_h::Int>());
  if (refine_opts.contains("Verbosity"))
    orSrv->SetVerbosity(refine_opts["Verbosity"].as<std::string>());
  if (refine_opts.contains("Length Histogram Min"))
    orSrv->SetLengthHistogramMin(
        refine_opts["Length Histogram Min"].as<Omega_h::Real>());
  if (refine_opts.contains("Length Histogram Max"))
    orSrv->SetLengthHistogramMax(
        refine_opts["Length Histogram Max"].as<Omega_h::Real>());
  if (refine_opts.contains("N length Histogram Bins"))
    orSrv->SetNlengthHistogramBins(
        refine_opts["N length Histogram Bins"].as<Omega_h::Int>());
  if (refine_opts.contains("N quality Histogram Bins"))
    orSrv->SetNqualityHistogramBins(
        refine_opts["N quality Histogram Bins"].as<Omega_h::Int>());
  if (refine_opts.contains("Should Refine"))
    orSrv->SetShouldRefine(refine_opts["Should Refine"].as<bool>());
  if (refine_opts.contains("Should Coarsen"))
    orSrv->SetShouldCoarsen(refine_opts["Should Coarsen"].as<bool>());
  if (refine_opts.contains("Should Swap"))
    orSrv->SetShouldSwap(refine_opts["Should Swap"].as<bool>());
  if (refine_opts.contains("Should Coarsen Slivers"))
    orSrv->SetShouldCoarsenSlivers(
        refine_opts["Should Coarsen Slivers"].as<bool>());
  if (refine_opts.contains("Should Prevent Coarsen Flip"))
    orSrv->SetShouldPreventCoarsenFlip(
        refine_opts["Should Prevent Coarsen Flip"].as<bool>());

  // Add transfer options
  if (refine_opts.contains("Transfer Options"))
    for (const auto &xfer : refine_opts["Transfer Options"].array_range())
      orSrv->AddTransferOpts(xfer["Name"].as<std::string>(),
                             Transfer(xfer["Method"].as<std::string>()),
                             xfer.get_with_default("Integral Name", ""));
  if (refine_opts.contains("Transfer Integral Options"))
    for (const auto &xfer :
         refine_opts["Transfer Integral Options"].array_range())
      orSrv->AddTransferOptsIntegralDiffuse(
          xfer["Integral Name"].as<std::string>(),
          xfer["Type"].as<std::string>(), xfer["Tolerance"].as<Omega_h::Real>(),
          xfer["Floor"].as<Omega_h::Real>());

  orSrv->Update();

  // Populate the output mesh
  out_mesh->takeGeoMesh(orSrv->GetOutput());
}

omegahRefineDriver *omegahRefineDriver::readJSON(
    const jsoncons::json &inputjson) {
  std::string ifName =
      inputjson["Mesh File Options"]["Input Mesh File"].as<std::string>();
  std::string ofName =
      inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  bool reconstructGeo = false;
  if (inputjson["Mesh File Options"].contains("Reconstruct Geometry"))
    reconstructGeo =
        inputjson["Mesh File Options"]["Reconstruct Geometry"].as<bool>();
  jsoncons::json refineOpts = inputjson["Refinement Options"];

  auto inMesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(ifName));
  if (reconstructGeo) {
    inMesh->reconstructGeo();
  }
  auto outMesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::New(ofName));
  auto *drv = new omegahRefineDriver(inMesh, outMesh, refineOpts);

  outMesh->write(ofName);

  return drv;
}

}  // namespace DRV
}  // namespace NEM
