#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "getDicts.H"

using namespace Foam;

std::unique_ptr<Foam::dictionary> getDicts::createControlDict(
    const bool write, const std::string &pth, const double &startT,
    const double &dt, const double &endT) {
  if (!Foam::isDir(pth) && (write)) Foam::mkDir(pth);
  std::unique_ptr<Foam::dictionary> controlDict_;
  controlDict_ =
      std::unique_ptr<Foam::dictionary>(new Foam::dictionary("controlDict"));
  Foam::dictionary fmfle("FoamFile");
  fmfle.add("version", "2.0");
  fmfle.add("format", "ascii");
  fmfle.add("class", "dictionary");
  fmfle.add("location", "\"system\"");
  fmfle.add("object", "controlDict");
  controlDict_->add("FoamFile", fmfle);
  controlDict_->add("startFrom", "startTime");
  controlDict_->add("startTime", startT);
  controlDict_->add("stopAt", "endTime");
  controlDict_->add("endTime", endT);
  controlDict_->add("deltaT", dt);
  controlDict_->add("writeControl", "timeStep");
  controlDict_->add("writeInterval", 1);
  controlDict_->add("purgeWrite", 0);
  controlDict_->add("writeFormat", "ascii");
  controlDict_->add("writePrecision", 6);
  controlDict_->add("writeCompression", "off");
  controlDict_->add("timeFormat", "general");
  controlDict_->add("timePrecision", 6);
  controlDict_->add("runTimeModifiable", "true");

  if (write) writeDict(pth + "controlDict", controlDict_);
  return controlDict_;
}

std::unique_ptr<Foam::dictionary> getDicts::createFvSchemes(
    const bool write, const std::string &pth) {
  if (!Foam::isDir(pth) && (write)) Foam::mkDir(pth);
  std::unique_ptr<Foam::dictionary> fvSchemes_;
  fvSchemes_ =
      std::unique_ptr<Foam::dictionary>(new Foam::dictionary("fvSchemes"));
  Foam::dictionary gradSchemes("gradSchemes");
  Foam::dictionary divSchemes("divSchemes");
  Foam::dictionary laplacianSchemes("laplacianSchemes");
  gradSchemes.add("default", "Gauss linear");
  gradSchemes.add("grad(p)", "Gauss linear");
  divSchemes.add("default", "none");
  divSchemes.add("div(phi,U)", "Gauss linear");
  laplacianSchemes.add("default", "none");
  laplacianSchemes.add("laplacian(nu,U)", "Gauss linear corrected");
  laplacianSchemes.add("laplacian((1|A(U)),p)", "Gauss linear corrected");
  Foam::dictionary fmfleFvscheme("FoamFile");
  fmfleFvscheme.add("version", "2.0");
  fmfleFvscheme.add("format", "ascii");
  fmfleFvscheme.add("class", "dictionary");
  fmfleFvscheme.add("location", "\"system\"");
  fmfleFvscheme.add("object", "fvSchemes");
  fvSchemes_->add("FoamFile", fmfleFvscheme);
  fvSchemes_->add("gradSchemes", gradSchemes);
  fvSchemes_->add("divSchemes", divSchemes);
  fvSchemes_->add("laplacianSchemes", laplacianSchemes);
  if (write) writeDict(pth + "fvSchemes", fvSchemes_);
  return fvSchemes_;
}

std::unique_ptr<Foam::dictionary> getDicts::createFvSolution(
    const bool write, const std::string &pth) {
  if (!Foam::isDir(pth) && (write)) Foam::mkDir(pth);
  std::unique_ptr<Foam::dictionary> fvSolution_;
  fvSolution_ =
      std::unique_ptr<Foam::dictionary>(new Foam::dictionary("fvSolution"));
  Foam::dictionary fmfleFvsol("FoamFile");
  fmfleFvsol.add("version", "2.0");
  fmfleFvsol.add("format", "ascii");
  fmfleFvsol.add("class", "dictionary");
  fmfleFvsol.add("location", "\"system\"");
  fmfleFvsol.add("object", "fvSolution");
  fvSolution_->add("FoamFile", fmfleFvsol);
  if (write) writeDict(pth + "fvSolution", fvSolution_);
  return fvSolution_;
}

std::unique_ptr<Foam::dictionary> getDicts::createDynamicMeshDict(
    const bool write, const std::string &pth) {
  if (!Foam::isDir(pth) && (write)) Foam::mkDir(pth);
  std::unique_ptr<Foam::dictionary> dynMshDict;
  dynMshDict = std::unique_ptr<Foam::dictionary>(
      new Foam::dictionary("dynamicMeshDict"));
  Foam::dictionary fmfleFvsol("FoamFile");
  fmfleFvsol.add("version", "2.0");
  fmfleFvsol.add("format", "ascii");
  fmfleFvsol.add("class", "dictionary");
  fmfleFvsol.add("location", "\"system\"");
  fmfleFvsol.add("object", "fvSolution");
  dynMshDict->add("FoamFile", fmfleFvsol);

  std::string contText = "dynamicFvMesh   dynamicRefineFvMesh;\n\n";
  contText = contText + "refineInterval  1;\n\n";
  contText = contText + "field  ";
  contText = contText + "meshField;\n\n";
  contText = contText + "lowerRefineLevel  0;\n\n";
  contText = contText + "upperRefineLevel  1;\n\n";
  contText = contText + "unrefineLevel  10;\n\n";

  contText = contText + "nBufferLayers   3;\n\n";
  contText = contText + "maxRefinement   3;\n\n";
  contText = contText + "maxCells        2000000;\n\n";
  contText = contText + "correctFluxes\n";
  contText = contText + "(\n";
  contText = contText + ");\n";
  contText = contText + "dumpLevel       true;";
  Foam::dictionary tmptmpDuc =
      new Foam::dictionary(Foam::IStringStream(contText)(), true);
  dynMshDict->merge(tmptmpDuc);
  if (write) writeDict(pth + "dynamicMeshDict", dynMshDict);
  return dynMshDict;
}

void getDicts::writeDict(const std::string &path,
                         std::unique_ptr<Foam::dictionary> &dictToWrite) {
  Foam::fileName fDict_ = path;
  if (!Foam::exists(fDict_)) {
    Foam::OFstream outDict_(fDict_);
    IOobject::writeBanner(outDict_);
    dictToWrite->write(outDict_, false);
  }
}

void getDicts::writeBasicDicts(const std::string &pth, const double &startT,
                               const double &dt, const double &endT) {
  if (!Foam::isDir(pth)) Foam::mkDir(pth);
  createControlDict(true, pth, startT, dt, endT);
  createFvSchemes(true, pth);
  createFvSolution(true, pth);
}