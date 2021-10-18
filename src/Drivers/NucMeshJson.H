#ifndef NEMOSYS_NUCMESHJSON_H_
#define NEMOSYS_NUCMESHJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/NucMeshDriver.H"

#include <map>
#include <memory>
#include <string>
#include <type_traits>

#include "NucMesh/CirclesAndPolys.H"
#include "NucMesh/HexagonalArray.H"
#include "NucMesh/PolarArray.H"
#include "NucMesh/RectangularArray.H"
#include "NucMesh/ShapeBase.H"
#include "NucMesh/ShapesArray.H"

namespace NEM {
namespace DRV {
namespace JSON {
// primary template.
template <typename T>
struct function_traits : function_traits<decltype(&T::operator())> {};

// function
template <typename R, typename... Args>
struct function_traits<R(Args...)> {
  using result_type = R;
  using argument_types = std::tuple<Args...>;
};

// function pointers
template <typename R, typename... Args>
struct function_traits<R (*)(Args...)> {
  using result_type = R;
  using argument_types = std::tuple<Args...>;
};

// member functions
template <typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...)> {
  using result_type = R;
  using argument_types = std::tuple<Args...>;
};

template <typename T>
using first_argument_type = typename std::tuple_element<
    0, typename function_traits<T>::argument_types>::type;

static constexpr auto shapeType = "Type";
static constexpr auto arrayPattern = "Pattern";
static constexpr auto arrayShapes = "Shapes";

template <typename T, typename R, typename... Args>
void setter_helper(T &obj, R (T::*setter)(Args...), const jsoncons::json &json,
                   jsoncons::string_view key, bool checkContains) {
  if (json.contains(key)) {
    (obj.*setter)(json.at(key)
                      .as<typename std::decay<
                          first_argument_type<decltype(setter)>>::type>());
  } else if (checkContains) {
    throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                  "Not a NEM::NUCMESH::ShapeBase");
  }
}

template <typename T>
std::unique_ptr<NEM::NUCMESH::ShapeBase> construct_helper(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes);

inline std::unique_ptr<NEM::NUCMESH::ShapeBase> as_shape(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes);

template <typename T>
typename std::enable_if<std::is_base_of<NEM::NUCMESH::ShapeBase, T>::value,
                        void>::type
modify_helper(
    const jsoncons::json &json, T &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::ShapeBase &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>> &,
    bool) {
  if (json.contains("Center")) {
    auto centerSize = json.at("Center").array_value().size();
    if (centerSize == 3) {
      shape.setCenter(json.at("Center").as<std::array<double, 3>>());
    } else if (centerSize == 5) {
      auto center = json.at("Center").as<std::array<double, 5>>();
      shape.setCenter(NEM::NUCMESH::ShapeBase::getRotatedPoint(
          {center[0], center[1], center[2]}, {center[3], center[4]}));
    }
  }
}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::Circles &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {
  setter_helper(shape, &NEM::NUCMESH::Circles::setRings, json, "Rings",
                checkRequired);
  modify_helper(json, dynamic_cast<NEM::NUCMESH::ShapeBase &>(shape),
                savedShapes, checkRequired);
}

template <>
std::unique_ptr<NEM::NUCMESH::ShapeBase>
construct_helper<NEM::NUCMESH::Circles>(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  auto shape = new NEM::NUCMESH::Circles{{}};
  std::unique_ptr<NEM::NUCMESH::ShapeBase> out{shape};
  modify_helper(json, *shape, savedShapes, true);
  return out;
}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::CirclesAndPolys &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {
  setter_helper(shape, &NEM::NUCMESH::CirclesAndPolys::setNumSides, json,
                "Number of Sides", checkRequired);
  setter_helper(shape, &NEM::NUCMESH::CirclesAndPolys::setRings, json, "Rings",
                checkRequired);
  modify_helper(json, dynamic_cast<NEM::NUCMESH::ShapeBase &>(shape),
                savedShapes, checkRequired);
}

template <>
std::unique_ptr<NEM::NUCMESH::ShapeBase>
construct_helper<NEM::NUCMESH::CirclesAndPolys>(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  auto shape = new NEM::NUCMESH::CirclesAndPolys{{}, {}};
  std::unique_ptr<NEM::NUCMESH::ShapeBase> out{shape};
  modify_helper(json, *shape, savedShapes, true);
  return out;
}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::ShapesArray &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {
  if (json.contains(arrayShapes) && json.at(arrayShapes).is_array()) {
    std::size_t idx = 0;
    for (auto &shapeJson : json.at(arrayShapes).array_range()) {
      shape.setPatternShape(idx, as_shape(shapeJson, savedShapes));
      ++idx;
    }
    for (; idx < shape.getNumPatternShapes(); ++idx) {
      shape.setPatternShape(idx, nullptr);
    }
  } else if (checkRequired) {
    throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                  "Not a NEM::NUCMESH::ShapesArray");
  }
  modify_helper(json, dynamic_cast<NEM::NUCMESH::ShapeBase &>(shape),
                savedShapes, checkRequired);
}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::PolarArray &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {
  setter_helper(shape, &NEM::NUCMESH::PolarArray::setRadius, json, "Radius",
                checkRequired);
  setter_helper(shape, &NEM::NUCMESH::PolarArray::setStartAngle, json,
                "Start Angle", checkRequired);
  setter_helper(shape, &NEM::NUCMESH::PolarArray::setEndAngle, json,
                "End Angle", checkRequired);
  setter_helper(shape, &NEM::NUCMESH::PolarArray::setRotateWithArray, json,
                "Rotate with Array", false);
  if (json.contains(arrayPattern)) {
    auto patternIdJson = json.at(arrayPattern).array_value();
    if (!patternIdJson.empty()) {
      auto patternIdIter = patternIdJson.begin();
      const auto patternIdEnd = patternIdJson.end();
      for (std::size_t patternIdx = 0; patternIdx < shape.getNumSubshapes();
           ++patternIdx) {
        shape.setPattern(patternIdx, patternIdIter->as<std::size_t>());
        ++patternIdIter;
        if (patternIdIter == patternIdEnd) {
          patternIdIter = patternIdJson.begin();
        }
      }
    }
  } else if (checkRequired) {
    throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                  "Not a NEM::NUCMESH::PolarArray");
  }
  modify_helper(json, dynamic_cast<NEM::NUCMESH::ShapesArray &>(shape),
                savedShapes, checkRequired);
}

template <>
std::unique_ptr<NEM::NUCMESH::ShapeBase>
construct_helper<NEM::NUCMESH::PolarArray>(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  auto shape = new NEM::NUCMESH::PolarArray{
      json.at("Pattern").array_value().size(), {}, {}, {}};
  std::unique_ptr<NEM::NUCMESH::ShapeBase> out{shape};
  modify_helper(json, *shape, savedShapes, true);
  return out;
}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::RectangularArray &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {
  setter_helper(shape, &NEM::NUCMESH::RectangularArray::setGridDistance, json,
                "Grid Distance", checkRequired);
  if (json.contains(arrayPattern)) {
    const auto gridSize = shape.getGridDims();
    auto y = gridSize[1] - 1;
    for (auto &patternRow : json.at(arrayPattern).array_range()) {
      std::size_t x = 0;
      for (auto &patternJson : patternRow.array_range()) {
        shape.setPattern(x, y, patternJson.as<std::size_t>());
        ++x;
      }
      --y;
    }
  } else if (checkRequired) {
    throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                  "Not a NEM::NUCMESH::RectangularArray");
  }
  modify_helper(json, dynamic_cast<NEM::NUCMESH::ShapesArray &>(shape),
                savedShapes, checkRequired);
}

template <>
std::unique_ptr<NEM::NUCMESH::ShapeBase>
construct_helper<NEM::NUCMESH::RectangularArray>(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  auto patternJson = json.at("Pattern").array_value();
  auto shape = new NEM::NUCMESH::RectangularArray{
      {patternJson.size(), patternJson[0].size()}, {}};
  std::unique_ptr<NEM::NUCMESH::ShapeBase> out{shape};
  modify_helper(json, *shape, savedShapes, true);
  return out;
}

template <>
void modify_helper(
    const jsoncons::json &json, NEM::NUCMESH::HexagonalArray &shape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes,
    bool checkRequired) {
  setter_helper(shape, &NEM::NUCMESH::HexagonalArray::setGridDistance, json,
                "Grid Distance", checkRequired);
  if (json.contains(arrayPattern)) {
    int x = 0;
    for (auto &patternRow : json.at(arrayPattern).array_range()) {
      int y = 0;
      for (auto &patternJson : patternRow.array_range()) {
        shape.setPatternRowCol(x, y, patternJson.as<std::size_t>());
        ++y;
      }
      ++x;
    }
  } else if (checkRequired) {
    throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                  "Not a NEM::NUCMESH::RectangularArray");
  }
  modify_helper(json, dynamic_cast<NEM::NUCMESH::ShapesArray &>(shape),
                savedShapes, checkRequired);
}

template <>
std::unique_ptr<NEM::NUCMESH::ShapeBase>
construct_helper<NEM::NUCMESH::HexagonalArray>(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  auto shape = new NEM::NUCMESH::HexagonalArray{
      json.at("Pattern").array_value()[0].size(), {}};
  std::unique_ptr<NEM::NUCMESH::ShapeBase> out{shape};
  modify_helper(json, *shape, savedShapes, true);
  return out;
}

template <typename T>
std::unique_ptr<NEM::NUCMESH::ShapeBase> construct_from_alias_helper(
    const jsoncons::json &json, const T *const aliasShape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  auto output = new T{*aliasShape};
  modify_helper(json, *output, savedShapes, false);
  return std::unique_ptr<NEM::NUCMESH::ShapeBase>{output};
}

inline std::unique_ptr<NEM::NUCMESH::ShapeBase> construct_from_alias(
    const jsoncons::json &json, const NEM::NUCMESH::ShapeBase *const aliasShape,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  if (auto circles = dynamic_cast<const NEM::NUCMESH::Circles *>(aliasShape)) {
    return construct_from_alias_helper(json, circles, savedShapes);
  } else if (auto circlesAndPolys =
                 dynamic_cast<const NEM::NUCMESH::CirclesAndPolys *>(
                     aliasShape)) {
    return construct_from_alias_helper(json, circlesAndPolys, savedShapes);
  } else if (auto polar =
                 dynamic_cast<const NEM::NUCMESH::PolarArray *>(aliasShape)) {
    return construct_from_alias_helper(json, polar, savedShapes);
  } else if (auto rect = dynamic_cast<const NEM::NUCMESH::RectangularArray *>(
                 aliasShape)) {
    return construct_from_alias_helper(json, rect, savedShapes);
  } else if (auto hex = dynamic_cast<const NEM::NUCMESH::HexagonalArray *>(
                 aliasShape)) {
    return construct_from_alias_helper(json, hex, savedShapes);
  }
  throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                "Not a NEM::NUCMESH::ShapeBase");
}

inline std::unique_ptr<NEM::NUCMESH::ShapeBase> as_shape(
    const jsoncons::json &json,
    const std::map<std::string, std::unique_ptr<NEM::NUCMESH::ShapeBase>>
        &savedShapes) {
  if (json.is_null()) {
    return nullptr;
  }
  static constexpr auto loadSavedStr = "Saved Object";
  if (json.contains(loadSavedStr)) {
    auto findIter = savedShapes.find(json.at(loadSavedStr).as_string());
    if (findIter != savedShapes.end() && findIter->second) {
      return construct_from_alias(json, findIter->second.get(), savedShapes);
    }
  }
  static constexpr auto typeStr = "Type";
  if (json.contains(typeStr)) {
    auto jsonType = json.at(typeStr).as_string_view();
    if (jsonType == "Circles") {
      return construct_helper<NEM::NUCMESH::Circles>(json, savedShapes);
    } else if (jsonType == "Circles And Polys") {
      return construct_helper<NEM::NUCMESH::CirclesAndPolys>(json, savedShapes);
    } else if (jsonType == "Polar Array") {
      return construct_helper<NEM::NUCMESH::PolarArray>(json, savedShapes);
    } else if (jsonType == "Rectangular Array") {
      return construct_helper<NEM::NUCMESH::RectangularArray>(json,
                                                              savedShapes);
    } else if (jsonType == "Hexagonal Array") {
      return construct_helper<NEM::NUCMESH::HexagonalArray>(json, savedShapes);
    }
  }
  throw jsoncons::convert_error(jsoncons::convert_errc::conversion_failed,
                                "Not a NEM::NUCMESH::ShapeBase");
}

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::NucMeshDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, "NucMesh Options"),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::NucMeshDriver::programType;
     }))

JSONCONS_ENUM_NAME_TRAITS(NEM::NUCMESH::RingMeshOption::MeshingType, (TRI, "T"),
                          (QUAD, "Q"), (STRUCT, "S"))

JSONCONS_ENUM_NAME_TRAITS(NEM::NUCMESH::PolyRing::ShapeType, (CIRCLE, "Circle"),
                          (POLYGON, "Poly"))

namespace jsoncons {
template <typename Json>
struct json_type_traits<Json, NEM::NUCMESH::RingMeshOption> {
  using value_type = NEM::NUCMESH::RingMeshOption;
  using allocator_type = typename Json::allocator_type;
  using char_type = typename Json::char_type;
  using string_view_type = typename Json::string_view_type;
  static constexpr auto meshingTypeStr = "Type";
  static constexpr auto elemsStr = "Number of Elems";
  static bool is(const Json &ajson) noexcept {
    if (!ajson.is_object()) return false;
    if (!ajson.contains(meshingTypeStr)) return false;
    return true;
  }
  static value_type as(const Json &ajson) {
    using namespace NEM::NUCMESH;
    using namespace NEM::DRV::JSON;
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a NEM::NUCMESH::RingMeshOption");
    auto meshingType =
        ajson.at(meshingTypeStr).template as<RingMeshOption::MeshingType>();
    switch (meshingType) {
      case RingMeshOption::MeshingType::TRI:
        return RingMeshOption::ApplyTriMesh();
      case RingMeshOption::MeshingType::QUAD:
        return RingMeshOption::ApplyQuadMesh();
      case RingMeshOption::MeshingType::STRUCT:
        if (ajson.contains(elemsStr)) {
          return RingMeshOption::ApplyStructuredMesh(
              ajson.at(elemsStr)
                  .template as<std::decay<first_argument_type<
                      decltype(RingMeshOption::ApplyStructuredMesh)>>::type>());
        }
    }
    throw convert_error(convert_errc::conversion_failed,
                        "Not a NEM::NUCMESH::RingMeshOption");
  }
};
template <>
struct is_json_type_traits_declared<NEM::NUCMESH::RingMeshOption>
    : public std::true_type {};

template <typename Json>
struct json_type_traits<Json, NEM::NUCMESH::Ring> {
  using value_type = NEM::NUCMESH::Ring;
  using allocator_type = typename Json::allocator_type;
  using char_type = typename Json::char_type;
  using string_view_type = typename Json::string_view_type;
  static constexpr auto radiusStr = "Radius";
  static constexpr auto meshTypeStr = "Mesh";
  static constexpr auto materialStr = "Material";
  static constexpr auto sidesetStr = "Sideset";
  static bool is(const Json &ajson) noexcept {
    return ajson.is_object() && ajson.contains(radiusStr) &&
           ajson.contains(meshTypeStr);
  }
  static value_type as(const Json &ajson) {
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a NEM::NUCMESH::Ring");
    value_type aval(
        ajson.at(radiusStr).template as<decltype(value_type::radius)>(),
        ajson.at(meshTypeStr).template as<decltype(value_type::meshType)>());
    if (ajson.contains(materialStr)) {
      jsoncons::json_traits_helper<Json>::set_udt_member(ajson, materialStr,
                                                         aval.material);
    }
    if (ajson.contains(sidesetStr)) {
      jsoncons::json_traits_helper<Json>::set_udt_member(ajson, sidesetStr,
                                                         aval.sideset);
    }
    return aval;
  }
};
template <>
struct is_json_type_traits_declared<NEM::NUCMESH::Ring>
    : public std::true_type {};

template <typename Json>
struct json_type_traits<Json, NEM::NUCMESH::PolyRing> {
  using value_type = NEM::NUCMESH::PolyRing;
  using allocator_type = typename Json::allocator_type;
  using char_type = typename Json::char_type;
  using string_view_type = typename Json::string_view_type;
  static constexpr auto radiusStr =
      json_type_traits<Json, NEM::NUCMESH::Ring>::radiusStr;
  static constexpr auto meshTypeStr =
      json_type_traits<Json, NEM::NUCMESH::Ring>::meshTypeStr;
  static constexpr auto materialStr =
      json_type_traits<Json, NEM::NUCMESH::Ring>::materialStr;
  static constexpr auto sidesetStr =
      json_type_traits<Json, NEM::NUCMESH::Ring>::sidesetStr;
  static constexpr auto shapeTypeStr = "Shape Type";
  static constexpr auto rotationStr = "Rotation";
  static bool is(const Json &ajson) noexcept {
    return ajson.template is<NEM::NUCMESH::Ring>() &&
           ajson.contains(meshTypeStr);
  }
  static value_type as(const Json &ajson) {
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a NEM::NUCMESH::PolyRing");
    value_type aval(
        ajson.at(shapeTypeStr).template as<decltype(value_type::shapeType)>(),
        ajson.at(radiusStr).template as<decltype(value_type::radius)>(),
        ajson.at(meshTypeStr).template as<decltype(value_type::meshType)>());
    if (ajson.contains(rotationStr)) {
      jsoncons::json_traits_helper<Json>::set_udt_member(ajson, rotationStr,
                                                         aval.rotation);
    }
    if (ajson.contains(materialStr)) {
      jsoncons::json_traits_helper<Json>::set_udt_member(ajson, materialStr,
                                                         aval.material);
    }
    if (ajson.contains(sidesetStr)) {
      jsoncons::json_traits_helper<Json>::set_udt_member(ajson, sidesetStr,
                                                         aval.sideset);
    }
    return aval;
  }
};
template <>
struct is_json_type_traits_declared<NEM::NUCMESH::PolyRing>
    : public std::true_type {};

template <>
struct json_type_traits<json, NEM::DRV::NucMeshDriver::Opts> {
  using Json = json;
  using value_type = NEM::DRV::NucMeshDriver::Opts;
  constexpr static auto savedObjectStr = "Saved Objects";
  constexpr static auto savedObjectNameStr = "Name";
  constexpr static auto extrudeStr = "Extrude";
  static bool is(const Json &ajson) noexcept {
    if (!ajson.contains(NEM::DRV::JSON::arrayShapes) ||
        !ajson.at(NEM::DRV::JSON::arrayShapes).is_array())
      return false;
    if (ajson.contains(savedObjectStr) && !ajson.at(savedObjectStr).is_array())
      return false;
    return true;
  }
  static value_type as(const Json &ajson) {
    using namespace NEM::DRV;
    using namespace NEM::NUCMESH;
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a NEM::DRV::NucMeshDriver::Opts");
    std::map<std::string, std::unique_ptr<ShapeBase>> savedObjects;
    if (ajson.contains(savedObjectStr)) {
      auto &savedObjectsJson = ajson.at(savedObjectStr);
      for (auto &savedObjectJson : savedObjectsJson.array_range()) {
        if (savedObjectJson.contains(savedObjectNameStr)) {
          savedObjects.emplace(
              savedObjectJson.at(savedObjectNameStr).as_string_view(),
              JSON::as_shape(savedObjectJson, savedObjects));
        } else {
          throw convert_error(convert_errc::conversion_failed,
                              "Not a NEM::DRV::NucMeshDriver");
        }
      }
    }
    value_type aval{};
    for (auto &shape : ajson.at(JSON::arrayShapes).array_range()) {
      aval.geometryAndMesh.emplace_back(JSON::as_shape(shape, savedObjects));
    }
    if (ajson.contains(extrudeStr)) {
      json_traits_helper<Json>::set_udt_member(ajson, extrudeStr,
                                               aval.extrudeSteps);
    }
    return aval;
  }
};
template <>
struct is_json_type_traits_declared<NEM::DRV::NucMeshDriver::Opts>
    : public std::true_type {};
}  // namespace jsoncons

#endif  // NEMOSYS_NUCMESHJSON_H_
