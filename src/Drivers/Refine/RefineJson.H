/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_REFINEJSON_H_
#define NEMOSYS_REFINEJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/Refine/OmegahRefineDriver.H"
#include "Drivers/Refine/RefineDriver.H"

#ifdef HAVE_CFMSH
#  define IF_CFMSH(a) a,
#  include "Drivers/Refine/FoamRefineDriver.H"
#else
#  define IF_CFMSH(a)
#endif

#if defined(MLAMR) && defined(HAVE_CFMSH)
#  define IF_MLAMR_CFMSH(a) a,
#else
#  define IF_MLAMR_CFMSH(a)
#endif

#if defined(HAVE_GMSH)
#  define IF_GMSH(a) a,
#  include "Drivers/Refine/SizeFieldRefineDriver.H"
#  include "Drivers/Refine/UniformRefineDriver.H"
#  include "Drivers/Refine/Z2RefineDriver.H"
#else
#  define IF_GMSH(a)
#endif

namespace NEM {
namespace DRV {
namespace JSON {

static constexpr auto refOpts = "Refinement Options";
static constexpr auto refOperator = "Refinement Operator";
static constexpr auto refMethod = "Refinement Method";
static constexpr auto refTransferData = "Transfer Data";

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

// RefineDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(
    NEM::DRV::RefineDriver,
    (IF_MLAMR_CFMSH(NEM::DRV::FoamMLRefineDriver)
     IF_CFMSH(NEM::DRV::FoamRefineDriver)
     IF_GMSH(NEM::DRV::SizeFieldRefineDriver)
     IF_GMSH(NEM::DRV::UniformRefineDriver)
     IF_GMSH(NEM::DRV::Z2RefineDriver)
     NEM::DRV::OmegahRefineDriver),
    NEM::DRV::NemDriver, 2, (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (, , NEM::DRV::JSON::refOpts, NEM_JSON_RDONLY_OVERRIDE,
     NEM_JSON_CHECK_KEY_ONLY),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::RefineDriver::programType;
     }))

#ifdef HAVE_CFMSH
#  ifdef MLAMR
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::FoamMLRefineDriver, NEM::DRV::RefineDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::refOpts, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::FoamMLRefineDriver::Opts, 2,
    (method, NEM::DRV::JSON::refMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::FoamMLRefineDriver::Opts::method;
     }),
    (refCriteria, NEM::DRV::JSON::refOperator, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::FoamMLRefineDriver::Opts::refCriteria;
     }),
    (refineInterval, "Refinement Interval"),
    (maxRefinement, "Maximum Refinement"), (nBufferLayers, "Buffer Layers"),
    (maxCells, "Max Cells"), (writeFieldData, "Write Field Data?"),
    (writeMesh, "Write Mesh?"), (writeRefHistory, "Write Refinement Data?"),
    (timeStep, "Time Step"), (endTime, "End Time"), (startTime, "Start Time"),
    (mlModel, "ML Kernel"))
#  endif

NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::FoamRefineDriver, NEM::DRV::RefineDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::refOpts, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::FoamRefineDriver::Opts, 3,
    (method, NEM::DRV::JSON::refMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::FoamRefineDriver::Opts::method;
     }),
    (refCriteria, NEM::DRV::JSON::refOperator),
    (inputFieldFile, "Input Field File"),
    (refineInterval, "Refinement Interval"),
    (maxRefinement, "Maximum Refinement"), (nBufferLayers, "Buffer Layers"),
    (maxCells, "Max Cells"), (writeFieldData, "Write Field Data?"),
    (writeMesh, "Write Mesh?"), (writeRefHistory, "Write Refinement Data?"),
    (timeStep, "Time Step"), (endTime, "End Time"), (startTime, "Start Time"),
    (lowerRefineLevel, "Lower Refinement Level"),
    (upperRefineLevel, "Upper Refinement Level"),
    (unrefineAbove, "Unrefinement Above"),
    (unrefineBelow, "Unrefinement Below"))

JSONCONS_ENUM_NAME_TRAITS(NEM::DRV::FoamRefineDriver::Opts::Criteria,
                          (VALUE, "Value"), (GRADIENT, "Gradient"))
#endif

// OmegahRefineDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::OmegahRefineDriver, NEM::DRV::RefineDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::refOpts, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::OmegahRefineDriver::Opts, 2,
    (method, NEM::DRV::JSON::refMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::OmegahRefineDriver::Opts::method;
     }),
    (MetricSources, "Metric Sources"), (reconstructGeo, "Reconstruct Geometry"),
    (Verbose, "Verbose"), (ShouldLimitLengths, "Should Limit Lengths"),
    (MaxLength, "Max Length"), (MinLength, "Min Length"),
    (ShouldLimitGradation, "Should Limit Gradation"),
    (MaxGradationRate, "Max Gradation Rate"),
    (GradationConvergenceTolerance, "Gradation Convergence Tolerance"),
    (ShouldLimitElementCount, "Should Limit Element Count"),
    (MaxElementCount, "Max ElementCount"),
    (MinElementCount, "Min Element Count"),
    (ElementCountOverRelaxation, "Element Count Over Relaxation"),
    (NsmoothingSteps, "N smoothing Steps"),
    (MinLengthDesired, "Min Length Desired"),
    (MaxLengthDesired, "Max Length Desired"),
    (MaxLengthAllowed, "Max Length Allowed"),
    (MinQualityAllowed, "Min Quality Allowed"),
    (MinQualityDesired, "Min Quality Desired"),
    (NsliverLayers, "N sliver Layers"), (Verbosity, "Verbosity"),
    (LengthHistogramMin, "Length Histogram Min"),
    (LengthHistogramMax, "Length Histogram Max"),
    (NlengthHistogramBins, "N length Histogram Bins"),
    (NqualityHistogramBins, "N quality Histogram Bins"),
    (ShouldRefine, "Should Refine"), (ShouldCoarsen, "Should Coarsen"),
    (ShouldSwap, "Should Swap"),
    (ShouldCoarsenSlivers, "Should Coarsen Slivers"),
    (ShouldPreventCoarsenFlip, "Should Prevent Coarsen Flip"),
    (TransferOpts, "Transfer Options"),
    (TransferOptsIntegralDiffuse, "Transfer Integral Options"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::SRV::omegahRefineMetricSource, 1,
                              (type, "Type"), (knob, "Knob"),
                              (tag_name, "Tag Name"), (isotropy, "Isotropy"),
                              (scales, "Scales"))

JSONCONS_ENUM_NAME_TRAITS(Omega_h_Source, (OMEGA_H_CONSTANT, "Constant"),
                          (OMEGA_H_VARIATION, "Variation"),
                          (OMEGA_H_DERIVATIVE, "Derivative"),
                          (OMEGA_H_GIVEN, "Given"),
                          (OMEGA_H_IMPLIED, "Implied"),
                          (OMEGA_H_CURVATURE, "Curvature"))

JSONCONS_ENUM_NAME_TRAITS(Omega_h_Isotropy,
                          (OMEGA_H_ANISOTROPIC, "Anisotropic"),
                          (OMEGA_H_ISO_LENGTH, "Length"),
                          (OMEGA_H_ISO_SIZE, "Size"))

JSONCONS_ENUM_NAME_TRAITS(Omega_h_Scales, (OMEGA_H_ABSOLUTE, "Absolute"),
                          (OMEGA_H_SCALES, "Scales"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::OmegahRefineDriver::Transfer, 2,
                              (arrayName, "Name"), (method, "Method"),
                              (integralName, "Integral Name"))

JSONCONS_ENUM_NAME_TRAITS(Omega_h_Transfer, (OMEGA_H_INHERIT, "Inherit"),
                          (OMEGA_H_LINEAR_INTERP, "Linear Interp"),
                          (OMEGA_H_METRIC, "Metric"),
                          (OMEGA_H_DENSITY, "Density"),
                          (OMEGA_H_CONSERVE, "Conserve"),
                          (OMEGA_H_MOMENTUM_VELOCITY, "Momentum Velocity"),
                          (OMEGA_H_POINTWISE, "Pointwise"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::OmegahRefineDriver::VarCompare, 1,
                              (integralName, "Integral Name"), (type, "Type"),
                              (tolerance, "Tolerance"), (floor, "Floor"))

// SizeFieldRefineDriver
#ifdef HAVE_GMSH
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::SizeFieldRefineDriver, NEM::DRV::RefineDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::refOpts, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::SizeFieldRefineDriver::Opts, 5,
                              (method, NEM::DRV::JSON::refMethod),
                              (arrayName, "Array Name"),
                              (stdDevMult, "StdDev Multiplier"),
                              (maxIsMin, "Max Is Min for Scaling"),
                              (transferData, NEM::DRV::JSON::refTransferData),
                              (sizeFactor, "Size Factor"))

JSONCONS_ENUM_NAME_TRAITS(NEM::DRV::SizeFieldRefineDriver::Opts::Method,
                          (VALUE, "value"), (GRADIENT, "gradient"))

// UniformRefineDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::UniformRefineDriver, NEM::DRV::RefineDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::refOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::UniformRefineDriver::Opts,
    (method, NEM::DRV::JSON::refMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::UniformRefineDriver::Opts::method;
     }),
    (transferData, NEM::DRV::JSON::refTransferData),
    (edgeScale, "Edge Scaling"))

// Z2RefineDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::Z2RefineDriver, NEM::DRV::RefineDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::refOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::Z2RefineDriver::Opts,
    (method, NEM::DRV::JSON::refMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::Z2RefineDriver::Opts::method;
     }),
    (transferData, NEM::DRV::JSON::refTransferData), (arrayName, "Array Name"),
    (order, "Order"))
#endif

#undef IF_CFMSH
#undef IF_MLAMR_CFMSH
#undef IF_GMSH

#endif  // NEMOSYS_REFINEJSON_H_
