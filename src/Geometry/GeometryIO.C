#include "Geometry/GeometryIO.H"

#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <IGESControl_Writer.hxx>
#include <STEPControl_Writer.hxx>
#include <StlAPI_Writer.hxx>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace GEO {

GeomType GeomTypeFromFilename(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);

  if (fileExt == ".step" || fileExt == ".stp") {
    return GeomType::STEP_GEO;
  } else if (fileExt == ".iges") {
    return GeomType::IGES_GEO;
  } else if (fileExt == ".stl") {
    return GeomType::STL_GEO;
  } else {
    std::cerr << "File extension " << fileExt << " is not supported."
              << std::endl;
    exit(1);
  }
}

GeometryIO::GeometryIO(GeoManager *geoManager, std::string &filename)
    : geoManager_(geoManager) {}

GeometryIO::~GeometryIO() = default;

TopoDS_Shape GeometryIO::sewShapes() {
  BRepBuilderAPI_Sewing sewer{};
  for (auto &shape : geoManager_->getMap()) {
    if (!geoManager_->isChild(shape.first)) {
      sewer.Add(shape.first);
    }
  }
  sewer.Perform();
  TopoDS_Shape sewnShape = sewer.SewedShape();
  return sewnShape;
}

void GeometryIO::write(std::string &filename) {
  GeomType geomType = GeomTypeFromFilename(filename);
  TopoDS_Shape sewnShape = sewShapes();
  if (geomType == GeomType::STEP_GEO)
    writeSTEP(filename, sewnShape);
  if (geomType == GeomType::IGES_GEO)
    writeIGES(filename, sewnShape);
  if (geomType == GeomType::STL_GEO)
    writeSTL(filename, sewnShape);
}

void GeometryIO::writeSTEP(std::string &filename, TopoDS_Shape &sewnShape) {
  // Create the writer object
  STEPControl_Writer writer;
  // Set the write mode
  STEPControl_StepModelType mode = STEPControl_AsIs;
  // STEPControl_StepModelType mode = STEPControl_ShellBasedSurfaceModel;
  //  Add the resultant sewed shape to the writer
  IFSelect_ReturnStatus status = writer.Transfer(sewnShape, mode);
  if (!status) {
    std::cerr << " - Error" << std::endl;
  }
  status = writer.Write(filename.c_str());
  if (status) {
    std::cout << " - File " << filename << " written successfully."
              << std::endl;
  } else {
    std::cerr << " - File " << filename << " writing failed." << std::endl;
  }
}

void GeometryIO::writeIGES(std::string &filename, TopoDS_Shape &sewnShape) {
  IGESControl_Writer writer;
  bool status = writer.AddShape(sewnShape);
  if (!status)
    std::cerr << "Error - shape not added to IGES writer" << std::endl;
  else {
    status = writer.Write(filename.c_str());
    if (status) {
      std::cout << " - File " << filename << " written successfully."
                << std::endl;
    } else {
      std::cerr << " - File " << filename << " writing failed." << std::endl;
    }
  }
}

void GeometryIO::writeSTL(std::string &filename, TopoDS_Shape &sewnShape) {
  // Need to call BRepMesh to triangulate the shapes before STL writing
  BRepMesh_IncrementalMesh(sewnShape, 0.1, Standard_True, 0.2, Standard_False);
  StlAPI_Writer writer;
  bool status = writer.Write(sewnShape, filename.c_str());
  if (status) {
    std::cout << " - File " << filename << " written successfully."
              << std::endl;
  } else {
    std::cout << " - File " << filename << " writing failed." << std::endl;
  }
}

/*
bool GeometryIO::isChild(const TopoDS_Shape &shape) const {
  int shapeDim;
  switch (shape.ShapeType()) {
    case TopAbs_COMPSOLID:
    case TopAbs_SOLID: shapeDim = 3; break;
    case TopAbs_SHELL:
    case TopAbs_FACE: shapeDim = 2; break;
    case TopAbs_WIRE:
    case TopAbs_EDGE: shapeDim = 1; break;
    case TopAbs_VERTEX: shapeDim = 0; break;
    case TopAbs_COMPOUND:
    case TopAbs_SHAPE:
    default: return false;
  }
  return shapeDim < geoManager_->getDim();
}
*/

}  // namespace GEO
}  // namespace NEM
