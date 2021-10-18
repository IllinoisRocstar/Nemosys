#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/inpGeoMesh.H"

#include <algorithm>
#include <fstream>
#include <jsoncons/json.hpp>
#include <jsoncons_ext/csv/csv.hpp>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "AuxiliaryFunctions.H"

#ifdef HAVE_GMSH
#  include <gmsh.h>
#endif

namespace NEM {
namespace MSH {

namespace {

constexpr std::array<std::tuple<VTKCellType, const char *, int>, 5> vtk2inp{
    {{VTK_TRIANGLE, "CPS3", 2},
     {VTK_QUAD, "CPS4", 2},
     {VTK_TETRA, "C3D4", 3},
     {VTK_HEXAHEDRON, "C3D8", 3},
     {VTK_WEDGE, "C3D6", 3}}};

int inpSide2vtkSide(VTKCellType cellType, int inpSide) {
  switch (cellType) {
    case VTK_TETRA:
      switch (inpSide) {
        case 1: return 3;
        case 2: return 0;
        case 3: return 1;
        case 4: return 2;
        default: return -1;
      }
    case VTK_HEXAHEDRON:
      switch (inpSide) {
        case 1: return 4;
        case 2: return 5;
        case 3: return 2;
        case 4: return 1;
        case 5: return 3;
        case 6: return 0;
        default: return -1;
      }
    default: return inpSide - 1;
  }
}

int vtkSide2inpSide(VTKCellType cellType, int vtkSide) {
  switch (cellType) {
    case VTK_TETRA:
      switch (vtkSide) {
        case 0: return 2;
        case 1: return 3;
        case 2: return 4;
        case 3: return 1;
        default: return -1;
      }
    case VTK_HEXAHEDRON:
      switch (vtkSide) {
        case 0: return 6;
        case 1: return 4;
        case 2: return 3;
        case 3: return 5;
        case 4: return 1;
        case 5: return 2;
        default: return -1;
      }
    default: return vtkSide + 1;
  }
}

class InpParser {
 public:
  InpParser() = default;

  /**
   * @brief Class to store the result of parser a .inp file
   */
  struct Mesh {
    struct Elem {
      Elem(vtkIdType id, VTKCellType cellType, std::vector<vtkIdType> points)
          : id(id), cellType(cellType), points(std::move(points)) {}
      Elem(vtkIdType id, VTKCellType cellType)
          : id(id), cellType(cellType), points() {}
      /**
       * @brief id in .inp file
       */
      vtkIdType id;
      VTKCellType cellType;
      /**
       * @brief points given by id in .inp file
       */
      std::vector<vtkIdType> points;
    };
    int maxDim{0};
    /**
     * @brief Each entry stores (node ID in .inp file, coordinates)
     */
    std::vector<std::pair<vtkIdType, std::vector<double>>> nodes{};
    /**
     * @brief If the ELEMENT keyword line contains ELSET=name, then the element
     * is stored in elems[name], otherwise, stored in elems[std::string{}].
     */
    std::map<std::string, std::vector<Elem>> elems{};
    /**
     * @brief Map from ELSET name to element ids (ids given by .inp file)
     */
    std::map<std::string, std::vector<vtkIdType>> elSets{};
    /**
     * @brief Map from NSET name to node ids (ids given by .inp file)
     */
    std::map<std::string, std::vector<vtkIdType>> nodeSets{};
    /**
     * @brief Map from SURFACE name to (element id, side) (id and side both use
     * .inp IDs)
     */
    std::map<std::string, std::vector<std::pair<vtkIdType, int>>> surfaces{};
  };

  /**
   * @brief Parse line (ended by newline in file)
   * @param line line of file
   */
  void parseLine(const std::string &line) {
    if (line.at(0) == '*') {
      if (!this->setState(parseKwLine(line))) {
        return;  // Do nothing if comment line
      };
    } else {
      if (this->currContinue == Continue::KEYWORD) {
        bool continueNextLine{false};
        auto parsed = parseCSV(line);
        for (const auto &field : parsed[0].array_range()) {
          if (field.is_null()) {
            continueNextLine = true;
            break;
          } else if (field.is_array()) {
            this->currKwParams.emplace(field[0].as_string(),
                                       field[1].as_string());
          } else {
            this->currKwParams.emplace(field.as_string(), std::string{});
          }
        }
        this->currContinue =
            continueNextLine ? Continue::KEYWORD : Continue::NONE;
      } else {
        // Pretend to parse sections w/ unknown keywords
        if (this->currKw == Keyword::UNKNOWN) { return; }
        jsoncons::json parsed = parseCSV(line);
        auto parsedRange = parsed.array_range();
        auto iter = parsedRange.begin();
        bool continueNextLine = (parsedRange.end() - 1)->is_null();
        if (this->currKw == Keyword::NODE) {
          if (this->currContinue != Continue::DATA) {
            this->currId = (iter++)->as<vtkIdType>();
            if (!this->currSet.empty()) {
              this->mesh_.nodeSets[currSet].emplace_back(this->currId);
            }
            this->mesh_.nodes.emplace_back();
            auto &back = this->mesh_.nodes.back();
            back.first = this->currId;
            back.second.reserve(3);
          }
          auto &back = this->mesh_.nodes.back();
          for (; iter != parsedRange.end(); ++iter) {
            if (!iter->is_null()) {
              back.second.emplace_back(iter->as_double());
            }
          }
        } else if (this->currKw == Keyword::ELEMENT) {
          if (this->currCellType != VTK_EMPTY_CELL) {
            auto &elSet = this->mesh_.elems[this->currSet];
            if (this->currContinue != Continue::DATA) {
              this->currId = (iter++)->as<vtkIdType>();
              elSet.emplace_back(this->currId, this->currCellType);
            }
            auto &back = elSet.back();
            for (; iter != parsedRange.end(); ++iter) {
              if (!iter->is_null()) {
                back.points.emplace_back(iter->as<vtkIdType>());
              }
            }
          }
        } else if (this->currKw == Keyword::SURFACE) {
          auto iterType = this->currKwParams.find("TYPE");
          if (iterType == this->currKwParams.end() ||
              iterType->second != "ELEMENT") {
            auto &surf = this->mesh_.surfaces[this->currSet];
            auto iterElSet = this->mesh_.elSets.find(iter->as_string());
            if (iterElSet != this->mesh_.elSets.end()) {
              ++iter;
              auto side = std::stoi(iter->as_string().substr(1));
              for (const auto &elem : iterElSet->second) {
                surf.emplace_back(elem, side);
              }
            } else {
              auto elem = (iter++)->as<vtkIdType>();
              auto side = std::stoi(iter->as_string().substr(1));
              surf.emplace_back(elem, side);
            }
          }
        } else if (this->currKw == Keyword::NSET ||
                   this->currKw == Keyword::ELSET) {
          auto iterGenerate = this->currKwParams.find("GENERATE");
          auto &set = this->currKw == Keyword::NSET
                          ? this->mesh_.nodeSets[this->currSet]
                          : this->mesh_.elSets[this->currSet];
          if (iterGenerate != this->currKwParams.end()) {
            auto lower = (iter++)->as<vtkIdType>();
            auto upper = (iter++)->as<vtkIdType>();
            auto stride = iter == parsedRange.end() ? 1 : (iter++)->as<int>();
            for (auto i = lower; i <= upper; i += stride) {
              set.emplace_back(i);
            }
          } else {
            for (; iter != parsedRange.end(); ++iter) {
              if (!iter->is_null()) {
                auto iterOtherSet =
                    this->currKw == Keyword::NSET
                        ? this->mesh_.nodeSets.find(iter->as_string())
                        : this->mesh_.elSets.find(iter->as_string());
                if (iterOtherSet != (this->currKw == Keyword::NSET
                                         ? this->mesh_.nodeSets.end()
                                         : this->mesh_.elSets.end())) {
                  set.insert(set.end(), iterOtherSet->second.begin(),
                             iterOtherSet->second.end());
                } else {
                  set.emplace_back(iter->as<vtkIdType>());
                }
              }
            }
          }
        }
        this->currContinue = continueNextLine ? Continue::DATA : Continue::NONE;
      }
    }
  };

 private:
  Mesh mesh_{};

 public:
  Mesh &getMesh() { return mesh_; }

 private:
  enum class Keyword {
    COMMENT,
    NODE,
    ELEMENT,
    ELSET,
    NSET,
    SURFACE,
    UNKNOWN,
  };

  // Current data, that has not been registered in the above data members
  // Includes data from the keyword line
  /**
   * @brief current keyword
   */
  Keyword currKw{Keyword::UNKNOWN};
  /**
   * @brief Params in current keyword line
   */
  std::map<std::string, std::string> currKwParams{};
  /**
   * @brief TYPE, if inside data lines following ELEMENT keyword
   */
  VTKCellType currCellType{VTK_EMPTY_CELL};
  enum class Continue { NONE, KEYWORD, DATA };
  /**
   * @brief If previous line was keyword line that is continued, KEYWORD. If
   * previous line was data line that is continued, DATA. Otherwose, NONE
   */
  Continue currContinue{Continue::NONE};
  /**
   * @brief If previous line was data line that is continued and current keyword
   * is NODE or ELEMEMT, the ID of the node/cell being processed. Otherwise,
   * undefined.
   */
  vtkIdType currId{-1};
  /**
   * @brief Name of ELSET/NSET/Surface (empty string for ELEMENT keyword that
   * does not specify ELSET)
   */
  std::string currSet{};

  /**
   * @brief Parse keyword line
   * @param line Line of file
   * @return (Keyword, All other parameters in keyword line, whether or not line
   * is continued on next line)
   */
  static std::tuple<Keyword, std::map<std::string, std::string>, bool>
  parseKwLine(const std::string &line) {
    static constexpr std::array<std::pair<Keyword, const char *>, 5> matchArr{
        {{Keyword::NODE, "*NODE"},
         {Keyword::ELEMENT, "*ELEMENT"},
         {Keyword::ELSET, "*ELSET"},
         {Keyword::NSET, "*NSET"},
         {Keyword::SURFACE, "*SURFACE"}}};
    if (line.rfind("**", 0) == 0) {
      return {Keyword::COMMENT, {}, false};
    } else {
      Keyword kw{Keyword::UNKNOWN};
      std::map<std::string, std::string> params{};
      auto parsed = parseCSV(line);
      auto range = parsed.array_range();
      auto iter = range.begin();
      auto keywordStr = iter->as_string();
      nemAux::toUpper(keywordStr);
      for (const auto &match : matchArr) {
        if (keywordStr == match.second) {
          kw = match.first;
          break;
        }
      }
      ++iter;
      bool continueNextLine = (range.end() - 1)->is_null();
      for (; iter != range.end(); ++iter) {
        if (!iter->is_null()) {
          if (iter->is_array()) {
            auto key = iter->at(0).as_string();
            nemAux::toUpper(key);
            params.emplace(std::move(key), iter->at(1).as_string());
          } else {
            auto key = iter->as_string();
            nemAux::toUpper(key);
            params[std::move(key)];
          }
        }
      }
      return {kw, std::move(params), continueNextLine};
    }
  };

  /**
   * @brief Sets the parser's state based on output of @c parseKwLine
   * @param kwLine output of @c parseKwLine
   * @return whether or not line is continued on next line
   */
  bool setState(
      std::tuple<Keyword, std::map<std::string, std::string>, bool> kwLine) {
    if (std::get<0>(kwLine) != Keyword::COMMENT) {
      this->currKw = std::get<0>(kwLine);
      this->currKwParams = std::move(std::get<1>(kwLine));
      this->currSet.clear();
      this->currCellType = VTK_EMPTY_CELL;
      if (this->currKw == Keyword::ELEMENT) {
        auto &cellTypeStr = this->currKwParams.at("TYPE");
        auto cellTypeIter = std::find_if(
            vtk2inp.begin(), vtk2inp.end(),
            [&cellTypeStr](const decltype(vtk2inp)::value_type &x) {
              auto cellType = std::get<1>(x);
              return std::equal(
                  cellType, cellType + std::char_traits<char>::length(cellType),
                  cellTypeStr.begin());
            });
        if (cellTypeIter != vtk2inp.end()) {
          auto dim = std::get<2>(*cellTypeIter);
          this->currCellType = std::get<0>(*cellTypeIter);
          if (this->mesh_.maxDim < dim) {
            this->mesh_.maxDim = dim;
            this->mesh_.elems.clear();
          } else if (this->mesh_.maxDim > dim) {
            // Ignore these elements
            this->currCellType = VTK_EMPTY_CELL;
          }
        } else {
          std::cerr << "Unsupported cell type " << cellTypeStr << '\n';
        }
        auto elSetParamIter = this->currKwParams.find("ELSET");
        if (elSetParamIter != this->currKwParams.end()) {
          this->currSet = elSetParamIter->second;
        }
      }
      if (this->currKw == Keyword::ELSET) {
        this->currSet = this->currKwParams.at("ELSET");
      }
      if (this->currKw == Keyword::NODE) {
        auto nsetIter = this->currKwParams.find("NSET");
        if (nsetIter != this->currKwParams.end()) {
          this->currSet = nsetIter->second;
        }
      }
      if (this->currKw == Keyword::NSET) {
        this->currSet = this->currKwParams.at("NSET");
      }
      if (this->currKw == Keyword::SURFACE) {
        this->currSet = this->currKwParams.at("NAME");
      }
      this->currContinue =
          std::get<2>(kwLine) ? Continue::KEYWORD : Continue::NONE;
      this->currId = -1;
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Extract comma-separated-values
   * @param line data line
   * @return JSON array - if line in file is continued, last entry will be null.
   * Keyword parameters split as arrays
   */
  static jsoncons::json parseCSV(const std::string &line) {
    return jsoncons::csv::decode_csv<jsoncons::json>(
        line, jsoncons::csv::basic_csv_options<char>{}
                  .assume_header(false)
                  .subfield_delimiter('=')
                  .unquoted_empty_value_is_null(true)
                  .trim(true))[0];
  }
};

template <typename T>
void writeCheckWidth(std::ostream &outStream, const T &data,
                     const std::string &delimiter = ", ",
                     std::size_t maxWidth = 80, std::size_t maxEntries = 16) {
  std::size_t width = 0;
  std::size_t entries = 0;
  auto lenDelimiter = delimiter.size();
  for (auto iter = data.begin(); iter != data.end(); ++iter) {
    bool last = iter == data.end() - 1;
    width += iter->size() + (last ? 0 : lenDelimiter);
    entries += 1;
    if (width > maxWidth || entries > maxEntries) {
      outStream << "\n";
      width = iter->size() + (last ? 0 : lenDelimiter);
      entries = 1;
    }
    outStream << *iter;
    if (last) {
      outStream << '\n';
    } else {
      outStream << ", ";
    }
  }
}

void writeCells(std::ostream &outStream, vtkDataSet *data,
                const std::string &geo, const std::string &groupArrName) {
  // ELEMENT
  std::map<std::pair<VTKCellType, int>, std::vector<vtkIdType>> elemBlocks{};
  auto entArr = data->GetCellData()->GetArray(groupArrName.c_str());
  for (vtkIdType i = 0; i < data->GetNumberOfCells(); ++i) {
    elemBlocks[{static_cast<VTKCellType>(data->GetCellType(i)),
                entArr ? static_cast<int>(entArr->GetComponent(i, 0)) : 0}]
        .emplace_back(i);
  }
#ifdef HAVE_GMSH
  if (!geo.empty()) {
    gmsh::model::setCurrent(geo);
  }
#endif
  for (const auto &elems : elemBlocks) {
    auto cellType = elems.first.first;
    auto inpType =
        std::find_if(vtk2inp.begin(), vtk2inp.end(),
                     [cellType](const decltype(vtk2inp)::value_type &x) {
                       return std::get<0>(x) == cellType;
                     });
    if (inpType == vtk2inp.end()) {
      std::cerr << "Unsupported element type" << cellType << '\n';
      continue;
    }
    std::string elSet{};
#ifdef HAVE_GMSH
    if (!geo.empty()) {
      gmsh::model::getEntityName(gmsh::model::getDimension(),
                                 elems.first.second, elSet);
    }
#endif
    outStream << "*ELEMENT, TYPE=" << std::get<1>(*inpType);
    if (!elSet.empty()) { outStream << " ELSET=" << elSet; }
    outStream << '\n';
    for (const auto &cellIdx : elems.second) {
      auto cell = data->GetCell(cellIdx);
      std::vector<std::string> outStr{};
      outStr.reserve(cell->GetNumberOfPoints() + 1);
      outStr.emplace_back(std::to_string(cellIdx + 1));
      if (cellType == VTK_WEDGE) {
        outStr.emplace_back(std::to_string(cell->GetPointId(0) + 1));
        outStr.emplace_back(std::to_string(cell->GetPointId(2) + 1));
        outStr.emplace_back(std::to_string(cell->GetPointId(1) + 1));
        outStr.emplace_back(std::to_string(cell->GetPointId(3) + 1));
        outStr.emplace_back(std::to_string(cell->GetPointId(5) + 1));
        outStr.emplace_back(std::to_string(cell->GetPointId(4) + 1));
      } else {
        for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
          outStr.emplace_back(std::to_string(cell->GetPointId(k) + 1));
        }
      }
      writeCheckWidth(outStream, outStr);
    }
  }
}

#ifdef HAVE_GMSH
/**
 * @brief Run @p func for each appearance of a cell in @p dataSet in each
 * physical group of the current gmsh model
 * @tparam Func Callable type like with signature void (const std::string &,
 * vtkIdType)
 * @param dataSet dataSet with cells
 * @param entArr array mapping cells in @p dataSet to gmsh entity
 * @param dim dimension of gmsh physical group
 * @param func Callable object signature like void (const std::string &,
 * vtkIdType)
 */
template <typename Func>
void eachGmshPhyGroup(vtkDataSet *dataSet, vtkDataArray *entArr, int dim,
                      Func &&func) {
  std::map<int, std::vector<std::string>> ent2PhysGroup;
  {
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, dim);
    for (const auto &dimTag : dimTags) {
      std::vector<int> phyGrpIds{};
      gmsh::model::getPhysicalGroupsForEntity(dimTag.first, dimTag.second,
                                              phyGrpIds);
      for (const auto &phyGrp : phyGrpIds) {
        auto &phyGrpNameVec = ent2PhysGroup[dimTag.second];
        phyGrpNameVec.emplace_back();
        gmsh::model::getPhysicalName(dimTag.first, phyGrp,
                                     phyGrpNameVec.back());
      }
    }
  }
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    for (const auto &phyGrp :
         ent2PhysGroup[static_cast<int>(entArr->GetComponent(i, 0))]) {
      func(phyGrp, i);
    }
  }
}
#endif

}  // namespace

vtkStandardNewMacro(inpGeoMesh)

inpGeoMesh *inpGeoMesh::Read(const std::string &fileName) {
  return new inpGeoMesh(fileName);
}

inpGeoMesh::inpGeoMesh(const std::string &fileName)
    : inpGeoMesh(inp2GM(fileName)) {}

inpGeoMesh::inpGeoMesh(std::pair<GeoMesh, InpSets> mesh)
    : geoMeshBase(std::move(mesh.first)), inpSets_(std::move(mesh.second)) {}

inpGeoMesh::inpGeoMesh() : inpGeoMesh(std::string{}) {}

void inpGeoMesh::write(const std::string &fileName) {
  const auto &gm = getGeoMesh();
  auto mesh = gm.mesh;
  std::ofstream outFile(fileName);
  {  // NODE
    outFile << "*NODE\n";
    auto points = mesh->GetPoints();
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
      double *point = points->GetPoint(i);
      std::array<std::string, 4> outStr{
          std::to_string(i + 1), std::to_string(point[0]),
          std::to_string(point[1]), std::to_string(point[2])};
      writeCheckWidth(outFile, outStr);
    }
  }
  // mesh (ELEMENT)
  writeCells(outFile, mesh, gm.geo, gm.link);
  // ELSET
  for (const auto &elSet : this->inpSets_.elSets) {
    if (!elSet.second.empty()) {
      outFile << "*ELSET, ELSET=" << elSet.first << '\n';
      std::vector<std::string> outStr;
      outStr.reserve(elSet.second.size());
      std::transform(elSet.second.begin(), elSet.second.end(),
                     std::back_inserter(outStr),
                     [](int x) { return std::to_string(x + 1); });
      writeCheckWidth(outFile, outStr);
    }
  }
  // NSET
  for (const auto &nodeSet : this->inpSets_.nodeSets) {
    if (!nodeSet.second.empty()) {
      outFile << "*NSET, NSET=" << nodeSet.first << '\n';
      std::vector<std::string> outStr;
      outStr.reserve(nodeSet.second.size());
      std::transform(nodeSet.second.begin(), nodeSet.second.end(),
                     std::back_inserter(outStr),
                     [](int x) { return std::to_string(x + 1); });
      writeCheckWidth(outFile, outStr);
    }
  }
  // SURFACE
  if (gm.sideSet.sides) {
    auto origCellArr = gm.sideSet.getOrigCellArr();
    auto cellFaceArr = gm.sideSet.getCellFaceArr();
    for (const auto &surface : this->inpSets_.surfaces) {
      if (!surface.second.empty()) {
        outFile << "*SURFACE, NAME=" << surface.first << ", TYPE=ELEMENT\n";
        for (const auto &id : surface.second) {
          auto origCellIdx = origCellArr->GetTypedComponent(id, 0);
          auto origCellType =
              static_cast<VTKCellType>(mesh->GetCellType(origCellIdx));
          outFile << origCellIdx + 1 << ", S"
                  << vtkSide2inpSide(origCellType,
                                     cellFaceArr->GetTypedComponent(id, 0))
                  << '\n';
        }
      }
    }
  }
}

void inpGeoMesh::report(std::ostream &out) const {
  this->geoMeshBase::report(out);
}

void inpGeoMesh::reconstructGeo() {
  this->geoMeshBase::reconstructGeo();
  this->resetNative();
}

void inpGeoMesh::resetNative() {
  this->inpSets_.nodeSets.clear();
  this->inpSets_.elSets.clear();
  this->inpSets_.surfaces.clear();
  // Default behavior:
  // Convert mesh entities to element sets and mesh physical groups to element
  // sets and node sets
  // Convert sideSet physical groups to surfaces and node sets
  auto gm = getGeoMesh();
#ifdef HAVE_GMSH
  if (!gm.geo.empty()) {
    gmsh::model::setCurrent(gm.geo);
    if (!gm.link.empty() && gm.mesh->GetNumberOfCells() > 0) {
      eachGmshPhyGroup(
          gm.mesh, gm.mesh->GetCellData()->GetArray(gm.link.c_str()),
          gm.mesh->GetCell(0)->GetCellDimension(),
          [this](const std::string &phyGroup, vtkIdType cellIdx) {
            this->inpSets_.elSets[phyGroup].emplace(cellIdx);
            auto cell = this->getGeoMesh().mesh->GetCell(cellIdx);
            for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
              this->inpSets_.nodeSets[phyGroup].emplace(cell->GetPointId(i));
            }
          });
    }
    if (gm.sideSet.sides && gm.sideSet.sides->GetNumberOfCells() > 0) {
      gm.findSide2OrigCell();
      eachGmshPhyGroup(
          gm.sideSet.sides, gm.sideSet.getGeoEntArr(),
          gm.sideSet.sides->GetCell(0)->GetCellDimension(),
          [this](const std::string &phyGroup, vtkIdType cellIdx) {
            this->inpSets_.surfaces[phyGroup].emplace(cellIdx);
            auto cell = this->getGeoMesh().sideSet.sides->GetCell(cellIdx);
            for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
              this->inpSets_.nodeSets[phyGroup].emplace(cell->GetPointId(i));
            }
          });
    }
  }
#endif
}

std::pair<geoMeshBase::GeoMesh, inpGeoMesh::InpSets> inpGeoMesh::inp2GM(
    const std::string &fileName) {
  auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  std::map<std::string, std::set<vtkIdType>> outNSets, outElSets, outSurfaces;
  SideSet sideSet{};
  std::string geoName{};
  std::string linkName{};
  if (fileName.empty()) { return {{mesh, {}, {}, {}}, {}}; }
  std::ifstream inFile(fileName);
  std::string line;
  InpParser parser{};
  while (std::getline(inFile, line)) { parser.parseLine(line); }
  auto parserMesh = parser.getMesh();

  std::map<vtkIdType, vtkIdType> nodeInpId2GMIdx;
  if (!parserMesh.nodes.empty()) {
    vtkNew<vtkPoints> points;
    points->Allocate(parserMesh.nodes.size());
    for (const auto &parserNode : parserMesh.nodes) {
      nodeInpId2GMIdx[parserNode.first] =
          points->InsertNextPoint(parserNode.second.data());
    }
    mesh->SetPoints(points);
  }
  auto nodeMapper = [&nodeInpId2GMIdx](vtkIdType x) {
    return nodeInpId2GMIdx.at(x);
  };
  for (const auto &nSet : parserMesh.nodeSets) {
    auto &outNSet = outNSets[nSet.first];
    std::transform(nSet.second.begin(), nSet.second.end(),
                   std::inserter(outNSet, outNSet.end()), nodeMapper);
  }
  std::map<vtkIdType, vtkIdType> inpId2GMIdx;
  if (!parserMesh.elems.empty()) {
    mesh->Allocate(parserMesh.elems.size());
#ifdef HAVE_GMSH
    GmshInterface::Initialize();
    geoName = "geoMesh_" + nemAux::getRandomString(6);
    gmsh::model::add(geoName);
    gmsh::model::setCurrent(geoName);
#endif
    linkName = GEO_ENT_DEFAULT_NAME;
    vtkNew<vtkIntArray> geoEntArr;
    geoEntArr->SetName(linkName.c_str());
#ifndef HAVE_GMSH
    int entTag = 1;
#endif
    for (const auto &elemSet : parserMesh.elems) {
#ifdef HAVE_GMSH
      auto entTag = gmsh::model::addDiscreteEntity(parserMesh.maxDim);
      gmsh::model::setEntityName(parserMesh.maxDim, entTag, elemSet.first);
#endif
      for (const auto &elem : elemSet.second) {
        auto points = elem.points;
        std::transform(points.begin(), points.end(), points.begin(),
                       nodeMapper);
        if (elem.cellType == VTK_WEDGE) {
          std::swap(points.at(1), points.at(2));
          std::swap(points.at(4), points.at(5));
        }
        inpId2GMIdx[elem.id] =
            mesh->InsertNextCell(elem.cellType, points.size(), points.data());
        geoEntArr->InsertNextValue(entTag);
      }
#ifndef HAVE_GMSH
      ++entTag;
#endif
    }
    mesh->GetCellData()->AddArray(geoEntArr);
  }
  for (const auto &elSet : parserMesh.elSets) {
    auto &outElSet = outElSets[elSet.first];
    for (const auto &elem : elSet.second) {
      auto iterCellIdxMap = inpId2GMIdx.find(elem);
      if (iterCellIdxMap != inpId2GMIdx.end()) {
        outElSet.emplace(iterCellIdxMap->second);
      }
    }
  }
  if (!parserMesh.surfaces.empty()) {
    auto sideSetPD = vtkSmartPointer<vtkPolyData>::New();
    sideSetPD->SetPoints(mesh->GetPoints());
    sideSetPD->Allocate();
    vtkNew<vtkIntArray> entArr;
    vtkNew<vtkIdTypeArray> origCellArr;
    origCellArr->SetNumberOfComponents(2);
    vtkNew<vtkIntArray> cellFaceArr;
    cellFaceArr->SetNumberOfComponents(2);
#ifndef HAVE_GMSH
    int entTag = 1;
#endif
    for (const auto &surf : parserMesh.surfaces) {
      auto &outSurface = outSurfaces[surf.first];
#ifdef HAVE_GMSH
      auto entTag = gmsh::model::addDiscreteEntity(parserMesh.maxDim - 1);
      auto phyGroupTag =
          gmsh::model::addPhysicalGroup(parserMesh.maxDim - 1, {entTag});
      gmsh::model::setPhysicalName(parserMesh.maxDim - 1, phyGroupTag,
                                   surf.first);
#endif
      for (const auto &side : surf.second) {
        auto iterCellIdxMap = inpId2GMIdx.find(side.first);
        if (iterCellIdxMap != inpId2GMIdx.end()) {
          auto origCell = mesh->GetCell(iterCellIdxMap->second);
          auto vtkSide = inpSide2vtkSide(
              static_cast<VTKCellType>(origCell->GetCellType()), side.second);
          auto sideCell = origCell->GetCellDimension() == 2
                              ? origCell->GetEdge(vtkSide)
                              : origCell->GetFace(vtkSide);
          auto sideIdx = sideSetPD->InsertNextCell(sideCell->GetCellType(),
                                                   sideCell->GetPointIds());
          entArr->InsertTypedComponent(sideIdx, 0, entTag);
          origCellArr->InsertTuple2(sideIdx, iterCellIdxMap->second, -1);
          cellFaceArr->InsertTuple2(sideIdx, vtkSide, -1);
          outSurface.emplace(sideIdx);
        }
      }
#ifndef HAVE_GMSH
      ++entTag;
#endif
    }
    sideSet = {sideSetPD, entArr, origCellArr, cellFaceArr};
  }
  return {{mesh, std::move(geoName), std::move(linkName), sideSet},
          {std::move(outElSets), std::move(outNSets), std::move(outSurfaces)}};
}

}  // namespace MSH
}  // namespace NEM
