%module pyNemosys

%include <stdint.i>  // For std::int32_t, std::int64_t, etc
// Unlike the other swig library headers, std_array.i does not #include <array>
%{
#include <array>
%}
%include <std_array.i>
%include <std_map.i>
%include <std_pair.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_set.i>
%include <std_shared_ptr.i>

%include "optional.i"

// SWIG cannot tell if a class has an implicitly declared copy ctor. All Driver classes should be copyable, so instruct
// SWIG to generate copy ctors in Python.
%copyctor;
// Instruct SWIG to wrap nested classes. See "Namespaces and Inner Classes" in the dev guide
%feature("flatnested");

// Rename (for SWIG) once-nested classes Foo::Bar in namespace BAZ (thus fully qualified name ::BAZ::Foo::Bar) to
// _Foo_Bar (trailing underscore to follow Python naming convention of private variables)
%define HIDE_NESTED_IN_NS(namespace)
%rename("%(regex:/^namespace::(.*)::(.*)$/_\\1_\\2/)s", %$isclass, %$isnested, fullname=1, regextarget=1) "^namespace::[^:]*::[^:]*$";
%enddef
HIDE_NESTED_IN_NS(NEM::DRV);
// Prepend nested classes outside any namespace with an underscore.
%rename("%(regex:/^(.*)::(.*)$/_\\1_\\2/)s", %$isclass, %$isnested, %$not %$innamespace, fullname=1, regextarget=1) "^[^:]*::[^:]*$";

// Variadic macro to rename nested types (innerClass and otherInnerClasses) of outerclass so that Foo::Bar (which were
// renamed _Foo_Bar above) is accessible as Foo.Bar in Python
%define EXPOSE_FLATNESTED_PY_CLASSES(outerclass, innerclass, otherInnerClasses...)
%pythoncode {
outerclass.innerclass = _ ## outerclass ## _ ## innerclass
};
#if #otherInnerClasses != ""
EXPOSE_FLATNESTED_PY_CLASSES(outerclass, otherInnerClasses);
#endif
%enddef

// Define these so that SWIG ignores them
#define NEMOSYS_EXPORT
#define NEMOSYS_NO_EXPORT
#define NEMOSYS_DEPRECATED_EXPORT
#define JSONCONS_TYPE_TRAITS_FRIEND
#define vtkTypeMacro(a, b)
#define vtkSetMacro(a, b)
#define vtkGetMacro(a, b)
#define vtkBooleanMacro(a, b)

%template(ArrArrDouble_3_3) std::array<std::array<double, 3>, 3>;
%template(ArrDouble_3) std::array<double, 3>;
%template(ArrInt_3) std::array<int, 3>;
%template(ArrStr_3) std::array<std::string, 3>;
%template(MapArrInt_3Str) std::map<std::array<int, 3>, std::string>;
%template(MapStrStr) std::map<std::string, std::string>;
%template(PairArrDouble_3) std::pair<std::array<double, 3>, std::array<double, 3>>;
%template(PairStrDouble) std::pair<std::string, double>;
%template(PairStrVecDouble) std::pair<std::string, std::vector<double>>;
%template(PairStrVecStr) std::pair<std::string, std::vector<std::string>>;
%template(SetStr) std::set<std::string>;
%template(VecDouble) std::vector<double>;
%template(VecInt) std::vector<int>;
%template(VecPairStrDouble) std::vector<std::pair<std::string, double>>;
%template(VecPairStrVecDouble) std::vector<std::pair<std::string, std::vector<double>>>;
%template(VecPairStrVecStr) std::vector<std::pair<std::string, std::vector<std::string>>>;
%template(VecStr) std::vector<std::string>;
%template(VecVecInt) std::vector<std::vector<int>>;

namespace jsoncons {
template <typename T>
class optional;
};

// From optional.i; See section "jsoncons::optional" in dev guide
OPTIONAL_TEMPLATE_IMMUTABLE(OptionalBool, jsoncons::optional, bool)
OPTIONAL_TEMPLATE_IMMUTABLE(OptionalDouble, jsoncons::optional, double)
OPTIONAL_TEMPLATE_IMMUTABLE(OptionalInt, jsoncons::optional, int)
OPTIONAL_TEMPLATE_IMMUTABLE(OptionalInt32, jsoncons::optional, std::int32_t)
OPTIONAL_TEMPLATE_IMMUTABLE(OptionalStr, jsoncons::optional, std::string)
OPTIONAL_TEMPLATE(OptionalVecStr, jsoncons::optional, std::vector<std::string>)

// SWIG turns C++ enums into constants (enums inside a class become attributes of the class) with scoped enums
// having the enum class name and an underscore as a prefix to the enum member name
// get_renamed_enum removes these attributes for a SWIG-wrapped C++ scoped enum and returns a Python enum.IntEnum
%pythonbegin {
import json as _json

def get_renamed_enum(enum_class_name, outer_class=None):
    import enum
    prefix = enum_class_name + '_'
    if outer_class is None:
        mangled_enumerators = {k : v for k, v in globals().items() if k.startswith(prefix)}
        for k in mangled_enumerators.keys():
            del globals()[k]
    else:
        mangled_enumerators = {k : getattr(outer_class, k) for k in dir(outer_class) if k.startswith(prefix)}
        for k in mangled_enumerators.keys():
            delattr(outer_class, k)
    renamed_enumerators = {k[len(prefix):] : v for k, v in mangled_enumerators.items()}
    return enum.IntEnum(enum_class_name, renamed_enumerators)
};

// Replace a scoped enum class (which is nested in another wrapped type, outer_class) which SWIG has turned into a
// series of constants into a Python enum
%define RENAME_INNER_ENUM(outer_class, enum_class)
%pythoncode {
outer_class.enum_class = get_renamed_enum(#enum_class, outer_class)
};
%enddef

// Replace a scoped enum class which SWIG has turned into a series of constants into a Python enum
%define RENAME_ENUM(enum_class)
%pythoncode {
enum_class = get_renamed_enum(#enum_class)
};
%enddef

// Implements workaround described in dev guide section "Containers as Data Members"
// module : imported name of library that SWIG generates
// fullclass : C++ type name (including namespace qualifiers)
// renamedclass : SWIG renamed class
// setter : name of a setter function to create (not meant for use, but must be unique within the class)
// member : name of data member (assumed to be same in Python and C++)
// Developer's note: The setter would ideally take by const reference instead of value (SWIG won't pass an rvalue to the
// function anyways). But if the setter's signature is changed to:
//   void(const decltype(fullclass::member) &)
// then the automatically generated getter returns by const ref.
// And if instead we extend the class with our own getter:
//   decltype(fullclass::member) & getter() {
//     return $self->member;
//   };
// it seems that SWIG adds `const &` when parsing `decltype(fullclass::member)`.
// Also note this is somewhat fragile because of it is assumed that (1) we can access SWIG's getter using member.fget
// and (2) SWIG names extended member functions by joining the class and the function name with an underscore.
%define CUSTOM_PROPERTY(module, fullclass, renamedclass, setter, member)
%extend fullclass {
    void setter(decltype(fullclass::member) other) {
      $self->member = std::move(other);
    };
    %pythoncode {
      member = property(member.fget, module.renamedclass ## _ ## setter)
    };
};
%enddef

// For the public data members {member, otherMembers...}, of the types namespace::class::innerclass, the following
// variadic macros allow setters to implicitly call constructors, so that, e.g., if the data member Foo::bar has type
// std::vector<int>, then "Foo.bar = [1, 2, 3]" will work in Python, as opposed to requiring
// "Foo.bar = VecInt([1, 2, 3])".
%define CUSTOM_PROPERTIES_INNER_NS(namespace, class, innerclass, member, otherMembers...)
CUSTOM_PROPERTY(_pyNemosys, namespace::class::innerclass, _ ## class ## _ ## innerclass ##, _Set_ ## member, member)
#if #otherMembers != ""
CUSTOM_PROPERTIES_INNER_NS(namespace, class, innerclass, otherMembers)
#endif
%enddef

%define CUSTOM_PROPERTIES(class, member, otherMembers...)
CUSTOM_PROPERTY(_pyNemosys, class, class, _Set_ ## member, member);
#if #otherMembers != ""
CUSTOM_PROPERTIES(class, otherMembers)
#endif
%enddef

%define CUSTOM_PROPERTIES_NS(namespace, class, member, otherMembers...)
CUSTOM_PROPERTY(_pyNemosys, namespace::class, class, _Set_ ## member, member);
#if #otherMembers != ""
CUSTOM_PROPERTIES_NS(namespace, class, otherMembers)
#endif
%enddef

%define CUSTOM_PROPERTIES_INNER(class, innerclass, member, otherMembers...)
CUSTOM_PROPERTY(_pyNemosys, class::innerclass, _ ## class ## _ ## innerclass, _Set_ ## member, member)
#if #otherMembers != ""
CUSTOM_PROPERTIES_INNER(class, innerclass, otherMembers)
#endif
%enddef

%{
#include <Drivers/NemDriver.H>
%}
%ignore NEM::DRV::NemDriver::readJSON;
%include <Drivers/NemDriver.H>

#ifdef HAVE_EPIC
%{
#include <Drivers/InputGenDriver.H>
%}
// See NucMesh wrapping for details about the following workarounds for Opts being jsoncons::json
%ignore NEM::DRV::InputGenDriver::getOpts;
%ignore NEM::DRV::InputGenDriver::setOpts(jsoncons::json);
%ignore NEM::DRV::InputGenDriver::InputGenDriver(std::string, jsoncons::json);
%include <Drivers/InputGenDriver.H>
%rename("%s") NEM::DRV::InputGenDriver::getOpts;
%pythonprepend NEM::DRV::InputGenDriver::InputGenDriver{
    opts = _json.dumps(opts)
};
%feature("shadow") NEM::DRV::InputGenDriver::setOpts {
def setOpts(self, opts):
    return $action(self, _json.dumps(opts))
};
%feature("shadow") NEM::DRV::InputGenDriver::getOpts {
def getOpts(self):
    return _json.loads($action(self))
};
%extend NEM::DRV::InputGenDriver {
    InputGenDriver(std::string service, const std::string &opts) {
      return new NEM::DRV::InputGenDriver(std::move(service), jsoncons::json::parse(opts));
    }
    std::string getOpts() const {
      return $self->getOpts().to_string();
    }
    void setOpts(const std::string &opts) {
      $self->setOpts(jsoncons::json::parse(opts));
    }
};
%pythoncode {
if hasattr(InputGenDriver.__init__, '__annotations__'):
    InputGenDriver.__init__.__annotations__.pop('opts', None)
};
#endif  // HAVE_EPIC

#ifdef HAVE_HDF5
%{
#include <Drivers/ProteusDriver.H>
%}
%include <Drivers/ProteusDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(ProteusDriver, Files, Opts)
#endif  // HAVE_HDF5

#ifdef HAVE_TEMPLATE_MESH
%{
#include <Drivers/TemplateMeshDriver.H>
%}
%ignore NEM::DRV::TemplateMeshDriver::Opts::templateParams;
%include <Drivers/TemplateMeshDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(TemplateMeshDriver, Opts)
%feature("shadow") NEM::DRV::TemplateMeshDriver::Opts::_getTemplateParams {
def _getTemplateParams(self):
    return _json.loads($action(self))
};
%feature("shadow") NEM::DRV::TemplateMeshDriver::Opts::_setTemplateParams {
def _setTemplateParams(self, templateParams):
    return $action(self, _json.dumps(templateParams))
};
%extend NEM::DRV::TemplateMeshDriver::Opts {
    std::string _getTemplateParams() const {
      return $self->templateParams.to_string();
    }
    void _setTemplateParams(const std::string &opts) {
      $self->templateParams = (jsoncons::json::parse(opts));
    }
};
%pythoncode {
TemplateMeshDriver.Files = DriverOutFile
TemplateMeshDriver.Opts.templateParams = property(TemplateMeshDriver.Opts._getTemplateParams, TemplateMeshDriver.Opts._setTemplateParams)
};
#endif  // HAVE_TEMPLATE_MESH

%{
#include <Drivers/AutoVerificationDriver.H>
%}
%include <Drivers/AutoVerificationDriver.H>
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, AutoVerificationDriver, Opts, arrayIds);
EXPOSE_FLATNESTED_PY_CLASSES(AutoVerificationDriver, Files, Opts);

%{
#include <Drivers/Conversion/ConversionDriver.H>
%}
%ignore NEM::DRV::ConversionDriver::genExo;
%include <Drivers/Conversion/ConversionDriver.H>

#ifdef HAVE_CFMSH
%{
#include <Drivers/Conversion/FoamToMshConversionDriver.H>
%}
%ignore NEM::DRV::FoamToMshConversionDriver::Opts;
%include <Drivers/Conversion/FoamToMshConversionDriver.H>
%pythoncode {
FoamToMshConversionDriver.Files = DriverOutFile
};

%{
#include <Drivers/Conversion/FoamToVtkConversionDriver.H>
%}
%ignore NEM::DRV::FoamToVtkConversionDriver::Opts;
%include <Drivers/Conversion/FoamToVtkConversionDriver.H>
%pythoncode {
FoamToVtkConversionDriver.Files = DriverOutFile
};

%{
#include <Drivers/Conversion/VtkToFoamConversionDriver.H>
%}
%ignore NEM::DRV::VtkToFoamConversionDriver::Opts;
%include <Drivers/Conversion/VtkToFoamConversionDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(VtkToFoamConversionDriver, Files)
#endif  // HAVE_CFMSH

%{
#include <Drivers/Conversion/GmshToExoConversionDriver.H>
%}
// See std_container.i and https://github.com/swig/swig/issues/1176
%ignore std::vector<NEM::DRV::GmshToExoConversionDriver::MeshData>::vector(size_type);
%ignore std::vector<NEM::DRV::GmshToExoConversionDriver::MeshData>::resize(size_type);
%template(VecMeshData) std::vector<NEM::DRV::GmshToExoConversionDriver::MeshData>;
%ignore std::vector<NEM::DRV::GmshToExoConversionDriver::PostProcTask>::vector(size_type);
%ignore std::vector<NEM::DRV::GmshToExoConversionDriver::PostProcTask>::resize(size_type);
%template(VecPostProcTask) std::vector<NEM::DRV::GmshToExoConversionDriver::PostProcTask>;
%include <Drivers/Conversion/GmshToExoConversionDriver.H>
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, GmshToExoConversionDriver, MeshData, sideSetNames, elmBlkNames);
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, GmshToExoConversionDriver, Opts, meshData, tasks);
EXPOSE_FLATNESTED_PY_CLASSES(GmshToExoConversionDriver, PostProcTask, MeshData, Opts)

%{
#include <Drivers/Conversion/GmshToVtkConversionDriver.H>
%}
%ignore NEM::DRV::GmshToVtkConversionDriver::Opts;
%include <Drivers/Conversion/GmshToVtkConversionDriver.H>
%pythoncode {
GmshToVtkConversionDriver.Files = DriverInOutFiles
};

%{
#include <Drivers/Conversion/ManipExoConversionDriver.H>
%}
%template(VecCombineBlocks) std::vector<NEM::DRV::ManipExoConversionDriver::CombineBlocks>;
%include <Drivers/Conversion/ManipExoConversionDriver.H>
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, ManipExoConversionDriver, Opts, combineBlocks)
EXPOSE_FLATNESTED_PY_CLASSES(ManipExoConversionDriver, Opts, CombineBlocks)
%pythoncode {
ManipExoConversionDriver.Files = DriverInOutFiles
};

%{
#include <Drivers/Conversion/SmartConversionDriver.H>
%}
%ignore NEM::DRV::SmartConversionDriver::Opts;
%include <Drivers/Conversion/SmartConversionDriver.H>
%pythoncode {
SmartConversionDriver.Files = DriverInOutFiles
};

%{
#include <Drivers/Conversion/VtkHexToTetConversionDriver.H>
%}
%ignore NEM::DRV::VtkHexToTetConversionDriver::Opts;
%include <Drivers/Conversion/VtkHexToTetConversionDriver.H>
%pythoncode {
VtkHexToTetConversionDriver.Files = DriverInOutFiles
};

%{
#include <Drivers/Conversion/VtkToCobaltConversionDriver.H>
%}
%ignore NEM::DRV::VtkToCobaltConversionDriver::Opts;
%include <Drivers/Conversion/VtkToCobaltConversionDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(VtkToCobaltConversionDriver, Files)

%{
#include <Drivers/Conversion/VtkToPatranConversionDriver.H>
%}
%shared_ptr(NEM::DRV::VtkToPatranConversionDriver::BoundaryCond);
%shared_ptr(NEM::DRV::VtkToPatranConversionDriver::FaceBC);
%shared_ptr(NEM::DRV::VtkToPatranConversionDriver::NodeBC);
%template(VecBoundaryCond) std::vector<std::shared_ptr<NEM::DRV::VtkToPatranConversionDriver::BoundaryCond>>;
%include <Drivers/Conversion/VtkToPatranConversionDriver.H>;
EXPOSE_FLATNESTED_PY_CLASSES(VtkToPatranConversionDriver, BoundaryCond, FaceBC, NodeBC, Opts)
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, VtkToPatranConversionDriver, Opts, bcInfo, nodePatchPreference);
%pythoncode {
VtkToPatranConversionDriver.Files = DriverInOutFiles
};

%{
#include <Drivers/Conversion/VtkToPntConversionDriver.H>
%}
%template(BlockMap) std::vector<PNTMesh::blockType>;
%include <Drivers/Conversion/VtkToPntConversionDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(VtkToPntConversionDriver, Opts)
%pythoncode {
VtkToPntConversionDriver.Files = DriverInOutFiles
};
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, VtkToPntConversionDriver, Opts, elemBlockMap);
%rename("$ignore", %$isfunction, %$isglobal, regextarget=1, fullname=1) "PNTMesh::.*";
%ignore PNTMesh::pntMesh;
%include <Mesh/pntMesh.H>
CUSTOM_PROPERTIES_NS(PNTMesh, blockType, elmIds, srfBCTag, srfBCEleRef, glbSrfId, adjBlkId, adjElmId, adjRefId, eConn)
RENAME_ENUM(elementType)
RENAME_ENUM(surfaceBCTag)

%{
#include <Drivers/MeshGen/MeshGenDriver.H>
%}
%include <Drivers/MeshGen/MeshGenDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(MeshGenDriver, MeshGenFiles)

#ifdef HAVE_CFMSH
%{
#include <Drivers/MeshGen/BlockMeshMeshGenDriver.H>
%}
%ignore NEM::DRV::BlockMeshMeshGenDriver::Opts;
%include <Drivers/MeshGen/BlockMeshMeshGenDriver.H>
%pythoncode {
BlockMeshMeshGenDriver.Files = DriverOutFile
};
%shared_ptr(bmShape);
%shared_ptr(bmBox);
%shared_ptr(bmSphere);
%shared_ptr(bmCylTaperedCone);
OPTIONAL_TEMPLATE(OptionalBmAutoGenBox, jsoncons::optional, bmAutoGenBox)
%include <MeshGeneration/blockMeshParams.H>
CUSTOM_PROPERTIES(bmAutoGenBox, offset)
CUSTOM_PROPERTIES(bmBox, init, len, smplGrading, coordsBox, autoGenerate)
CUSTOM_PROPERTIES(bmSphere, center, sphrGrading)
CUSTOM_PROPERTIES(bmCylTaperedCone, centerCyl, cylGrading)
CUSTOM_PROPERTIES(blockMeshParams, nCells)

%{
#include <Drivers/MeshGen/CFMeshMeshGenDriver.H>
%}
%ignore NEM::DRV::CFMeshMeshGenDriver::Opts;
%include <Drivers/MeshGen/CFMeshMeshGenDriver.H>
%pythoncode {
CFMeshMeshGenDriver.Files = MeshGenDriver.MeshGenFiles
};
%template(VecCfmObjRef) std::vector<cfmObjRef>;
%template(VecCfmLclRefPatch) std::vector<cfmLclRefPatch>;
%template(VecCfmPtchBndLyr) std::vector<cfmPtchBndLyr>;
%template(VecCfmNewPatch) std::vector<cfmNewPatch>;
OPTIONAL_TEMPLATE(OptionalCfmBoundaryLayer, jsoncons::optional, cfmBoundaryLayer)
OPTIONAL_TEMPLATE(OptionalCfmSrfFeatEdge, jsoncons::optional, cfmSrfFeatEdge)
OPTIONAL_TEMPLATE(OptionalCfmMeshQual, jsoncons::optional, cfmMeshQual)
OPTIONAL_TEMPLATE(OptionalCfmRenBndry, jsoncons::optional, cfmRenBndry)
%include <MeshGeneration/cfmeshParams.H>
CUSTOM_PROPERTIES(cfmObjRef, params)
CUSTOM_PROPERTIES(cfmRenBndry, newPatches)
CUSTOM_PROPERTIES(cfmBoundaryLayer, blPatches)
CUSTOM_PROPERTIES(cfmeshParams, boundaryLayers, srfEdge, objRefLst, improveMeshQuality, refPatches, renBndry)

%{
#include <Drivers/MeshGen/SnappyMeshMeshGenDriver.H>
%}
%ignore NEM::DRV::SnappyMeshMeshGenDriver::Opts;
%include <Drivers/MeshGen/SnappyMeshMeshGenDriver.H>
%pythoncode {
SnappyMeshMeshGenDriver.Files = MeshGenDriver.MeshGenFiles
};
%shared_ptr(shmSearchableShape);
%shared_ptr(shmSearchableBox);
%shared_ptr(shmSearchableCylinder);
%shared_ptr(shmSearchableSphere);
%template(VecShmSearchableShape) std::vector<std::shared_ptr<shmSearchableShape>>;
%template(VecShmSTLDefinition) std::vector<shmSTLDefinition>;
%template(VecShmFeatureEdgeRef) std::vector<shmFeatureEdgeRef>;
%template(VecShmSurfRefine) std::vector<shmSurfRefine>;
%template(VecShmRegionRefine) std::vector<shmRegionRefine>;
%template(VecShmLayers) std::vector<shmLayers>;
%include <MeshGeneration/snappymeshParams.H>
CUSTOM_PROPERTIES(shmSearchableBox, minBound, maxBound)
CUSTOM_PROPERTIES(shmSearchableCylinder, axisPoint1, axisPoint2)
CUSTOM_PROPERTIES(shmSearchableSphere, center)
CUSTOM_PROPERTIES(shmGeo, stlPatchDefs, srchShape)
CUSTOM_PROPERTIES(shmCastMeshControls, locMesh, ftrEdge, surfRefs, geomRefs)
CUSTOM_PROPERTIES(shmLayerControls, layerVec)
#endif  // HAVE_CFMSH

#ifdef HAVE_NGEN
%{
#include <Drivers/MeshGen/NetgenMeshGenDriver.H>
%}
%ignore NEM::DRV::NetgenMeshGenDriver::Opts;
%include <Drivers/MeshGen/NetgenMeshGenDriver.H>
%pythoncode {
NetgenMeshGenDriver.Files = MeshGenDriver.MeshGenFiles
};
%include <MeshGeneration/netgenParams.H>
#endif  // HAVE_NGEN

%{
#include <Drivers/MeshGen/GmshMeshGenDriver.H>
%}
%ignore NEM::DRV::GmshMeshGenDriver::Opts;
%include <Drivers/MeshGen/GmshMeshGenDriver.H>
%pythoncode {
GmshMeshGenDriver.Files = MeshGenDriver.MeshGenFiles
};
%template(VecVolSizeField) std::vector<NEM::GEN::volSizeField>;
%template(SetTransfiniteBlock) std::set<NEM::GEN::TransfiniteBlock>;
%ignore NEM::GEN::gmshParams::getMeshExtensions;
%include <MeshGeneration/gmshParams.H>
CUSTOM_PROPERTIES_NS(NEM::GEN, volSizeField, params, num_list_params, strg_list_params)
CUSTOM_PROPERTIES_NS(NEM::GEN, TransfiniteBlock, axis, vert, type, coef)
CUSTOM_PROPERTIES_NS(NEM::GEN, gmshParams, sizeFields, color2groupMap, transfiniteBlocks)

%{
#include <Drivers/MeshQualityDriver.H>
%}
%ignore NEM::DRV::CheckMeshQualDriver::Opts;
#ifdef HAVE_CFMSH
%template(VecCfmeshQualityParams) std::vector<cfmeshQualityParams>;
%ignore NEM::DRV::OptimizeMeshQualDriver::Opts;
#endif
%include <Drivers/MeshQualityDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(CheckMeshQualDriver, Files)
#ifdef HAVE_CFMSH
%include <MeshQuality/cfmeshQualityParams.H>
#endif

%{
#include <Drivers/NucMeshDriver.H>
%}
%ignore NEM::DRV::NucMeshRunner;
// Instead of wrapping jsoncons::json, pass C++ <-> Python via strings by:
// (1) don't wrap the original functions (including the ctor)
// (2) instruct SWIG to extend the C++ class with getters and setters that take strings
// (3) instruct SWIG to inject some code around/before via shadow/pythonprepend that uses Python's json to convert
//     from Python data structures to strings.
%ignore NEM::DRV::NucMeshDriver::getGeometryAndMesh;
%ignore NEM::DRV::NucMeshDriver::setGeometryAndMesh(jsoncons::json);
%ignore NEM::DRV::NucMeshDriver::NucMeshDriver(NEM::DRV::NucMeshDriver::Files, jsoncons::json);
%include <Drivers/NucMeshDriver.H>
// The previous liness ignore the actual C++ member functions, but we don't want to ignore the below members with the
// same names
%rename("%s") NEM::DRV::NucMeshDriver::getGeometryAndMesh;
%pythonprepend NEM::DRV::NucMeshDriver::NucMeshDriver{
    geometryAndMesh = _json.dumps(geometryAndMesh)
};
%feature("shadow") NEM::DRV::NucMeshDriver::setGeometryAndMesh {
def setGeometryAndMesh(self, geometryAndMesh):
    return $action(self, _json.dumps(geometryAndMesh))
};
%feature("shadow") NEM::DRV::NucMeshDriver::getGeometryAndMesh {
def getGeometryAndMesh(self):
    return _json.loads($action(self))
};
%extend NEM::DRV::NucMeshDriver {
    NucMeshDriver(Files files, const std::string &geometryAndMesh) {
      return new NEM::DRV::NucMeshDriver(std::move(files), jsoncons::json::parse(geometryAndMesh));
    }
    std::string getGeometryAndMesh() const {
      return $self->getGeometryAndMesh().to_string();
    }
    void setGeometryAndMesh(const std::string &geometryAndMesh) {
      $self->setGeometryAndMesh(jsoncons::json::parse(geometryAndMesh));
    }
};
// SWIG does its best to add Python type annotations, but this will now be wrong, so just remove it.
%pythoncode {
if hasattr(NucMeshDriver.__init__, '__annotations__'):
    NucMeshDriver.__init__.__annotations__.pop('geometryAndMesh', None)
};

%{
#include <Drivers/PackMesh/PackMeshDriver.H>
%}
%include <Drivers/PackMesh/PackMeshDriver.H>

#ifdef HAVE_CFMSH
%{
#include <Drivers/PackMesh/HexPackMeshDriver.H>
%}
%include <Drivers/PackMesh/HexPackMeshDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(HexPackMeshDriver, Files, Opts)
%include <MeshManipulationFoam/MeshManipulationFoamParams.H>
EXPOSE_FLATNESTED_PY_CLASSES(MeshManipulationFoamParams, SurfaceLambdaMuSmooth, SplitMeshRegions, MergeMeshes, CreatePatch, FoamToSurface, SurfaceSplitByManifold)
CUSTOM_PROPERTIES_INNER(MeshManipulationFoamParams, SurfaceSplitByManifold, pckRegionNames)
#endif  // HAVE_CFMSH

%{
#include <Drivers/PackMesh/SurfacePackMeshDriver.H>
%}
OPTIONAL_TEMPLATE(OptionalCustomDomain, jsoncons::optional, CustomDomain)
OPTIONAL_TEMPLATE(OptionalPeriod3DOpts, jsoncons::optional, Periodic3DOpts)
%include <Drivers/PackMesh/SurfacePackMeshDriver.H>
RENAME_INNER_ENUM(SurfacePackMeshDriver, PhysGrpOpts)
EXPOSE_FLATNESTED_PY_CLASSES(SurfacePackMeshDriver, Files, CustomDomain, Periodic3DOpts, Opts)
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, SurfacePackMeshDriver, CustomDomain, initial, length)
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, SurfacePackMeshDriver, Periodic3DOpts, customDomain, transferMesh)
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, SurfacePackMeshDriver, Opts, periodic3DOpts)

%{
#include <Drivers/Refine/RefineDriver.H>
%}
%include <Drivers/Refine/RefineDriver.H>
%pythoncode {
RefineDriver.Files = DriverInOutFiles
};

#ifdef HAVE_CFMSH
%{
#include <Drivers/Refine/FoamRefineDriver.H>
%}
%include <Drivers/Refine/FoamRefineDriver.H>
#ifdef MLAMR
EXPOSE_FLATNESTED_PY_CLASSES(FoamMLRefineDriver, Opts)
#endif  // MLAMR
EXPOSE_FLATNESTED_PY_CLASSES(FoamRefineDriver, Opts)
RENAME_INNER_ENUM(FoamRefineDriver.Opts, Criteria)
#endif  // HAVE_CFMSH

%{
#include <Drivers/Refine/OmegahRefineDriver.H>
%}
%ignore NEM::SRV::omegahRefineSrv;
%include <Services/omegahRefineSrv.H>
// Copy-pasted to avoid other defs
// Note using enum class as opposed to enum to use the RENAME_ENUM macro;
// it seems that typemaps are created to/from Python int so no issues
// associating these values with the actual enum.
enum class Omega_h_Transfer {
  OMEGA_H_INHERIT,
  OMEGA_H_LINEAR_INTERP,
  OMEGA_H_METRIC,
  OMEGA_H_DENSITY,
  OMEGA_H_CONSERVE,
  OMEGA_H_MOMENTUM_VELOCITY,
  OMEGA_H_POINTWISE,
};
RENAME_ENUM(Omega_h_Transfer)
enum class Omega_h_Source {
  OMEGA_H_CONSTANT,
  OMEGA_H_VARIATION,
  OMEGA_H_DERIVATIVE,
  OMEGA_H_GIVEN,
  OMEGA_H_IMPLIED,
  OMEGA_H_CURVATURE,
};
RENAME_ENUM(Omega_h_Source)
enum class Omega_h_Isotropy {
  OMEGA_H_ANISOTROPIC,
  OMEGA_H_ISO_LENGTH,
  OMEGA_H_ISO_SIZE,
};
RENAME_ENUM(Omega_h_Isotropy)
enum class Omega_h_Scales {
  OMEGA_H_ABSOLUTE,
  OMEGA_H_SCALES,
};
RENAME_ENUM(Omega_h_Scales)
%template(VecOmegahRefineMetricSource) std::vector<NEM::SRV::omegahRefineMetricSource>;
%ignore std::vector<NEM::DRV::OmegahRefineDriver::Transfer>::vector(size_type);
%ignore std::vector<NEM::DRV::OmegahRefineDriver::Transfer>::resize(size_type);
%template(VecTransfer) std::vector<NEM::DRV::OmegahRefineDriver::Transfer>;
%ignore std::vector<NEM::DRV::OmegahRefineDriver::VarCompare>::vector(size_type);
%ignore std::vector<NEM::DRV::OmegahRefineDriver::VarCompare>::resize(size_type);
%template(VecOmegahRefineVarCompare) std::vector<NEM::DRV::OmegahRefineDriver::VarCompare>;
%include <Drivers/Refine/OmegahRefineDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(OmegahRefineDriver, Opts, Transfer, VarCompare)
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, OmegahRefineDriver, Opts, MetricSources, TransferOpts, TransferOptsIntegralDiffuse)

%{
#include <Drivers/Refine/SizeFieldRefineDriver.H>
%}
%include <Drivers/Refine/SizeFieldRefineDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(SizeFieldRefineDriver, Opts)
RENAME_INNER_ENUM(SizeFieldRefineDriver.Opts, Method)

%{
#include <Drivers/Refine/UniformRefineDriver.H>
%}
%include <Drivers/Refine/UniformRefineDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(UniformRefineDriver, Opts)

%{
#include <Drivers/Refine/Z2RefineDriver.H>
%}
%include <Drivers/Refine/Z2RefineDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(Z2RefineDriver, Opts)

%{
#include <Drivers/TransferDriver.H>
%}
%include <Drivers/TransferDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(TransferDriver, Files, Opts)
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, TransferDriver, Opts, arrayNames)

%{
#include <Mesh/meshBase.H>
%}
%ignore meshBase::CreateShared;
%ignore meshBase::CreateUnique;
%newobject meshBase::Create;
%ignore sortNemId_tVec_compare;
%include <Mesh/meshBase.H>

%{
#include <Transfer/TransferBase.H>
%}
%shared_ptr(TransferBase);
%include <Transfer/TransferBase.H>

// Cleanup to prevent exposing symbols
%pythoncode {
del get_renamed_enum
};
