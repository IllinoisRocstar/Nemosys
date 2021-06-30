#include "Drivers/Conversion/GmshToExoConversionDriver.H"

#include <vtkCell.h>
#include "AuxiliaryFunctions.H"
#include "exoMesh.H"
#include "meshSrch.H"

// helper functions

namespace NEM {
namespace DRV {

GmshToExoConversionDriver::PostProcTask::PostProcTask(std::string taskFile)
    : taskFile(std::move(taskFile)) {}

GmshToExoConversionDriver::MeshData::MeshData(std::string meshFile)
    : meshFileName(std::move(meshFile)) {}

GmshToExoConversionDriver::Opts::Opts(int numMeshes,
                                      std::vector<MeshData> meshData,
                                      bool needsPostProc)
    : numMeshes(numMeshes),
      meshData(std::move(meshData)),
      needsPostProc(needsPostProc) {}

GmshToExoConversionDriver::GmshToExoConversionDriver(Files file, Opts opts)
    : file_(std::move(file)), opts_(std::move(opts)) {}

GmshToExoConversionDriver::GmshToExoConversionDriver()
    : GmshToExoConversionDriver(Files{{}}, Opts{{}, {}, {}}) {}

const GmshToExoConversionDriver::Files &GmshToExoConversionDriver::getFiles()
    const {
  return file_;
}

void GmshToExoConversionDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

const GmshToExoConversionDriver::Opts &GmshToExoConversionDriver::getOpts()
    const {
  return opts_;
}

void GmshToExoConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void GmshToExoConversionDriver::execute() const {
  std::cout << "Converting to EXODUS II..." << std::endl;

  // sanity check
  if (this->opts_.numMeshes == 0) {
    std::cerr << "Error: At least one mesh should be provided!\n";
    exit(-1);
  }

  // starting conversion operation
  auto em = new NEM::MSH::EXOMesh::exoMesh(this->file_.outputFile);

  // reading meshes
  int ieb = 0;          // Element Block counter.
  int ins = 0;          // Node Set counter.
  int iss = 0;          // Side Set counter.
  int ndeIdOffset = 0;  // for element blocks and node sets
  int elmIdOffset = 0;  // for side sets
  for (int iMsh = 0; iMsh < this->opts_.numMeshes; iMsh++) {
    const auto &meshData = this->opts_.meshData[iMsh];
    // reading input mesh
    std::string fileExt = nemAux::find_ext(meshData.meshFileName);
    if (fileExt == ".g" || fileExt == ".e" || fileExt == ".exo" ||
        fileExt == ".exodus") {
      // if already exodus format, stitch and continue.
      NEM::MSH::EXOMesh::exoMesh em2;
      em2.read(meshData.meshFileName);

      // Allow changing element block names using key-value.
      for (const auto &kv : meshData.elmBlkNames)
        em2.setBlockName(std::stoi(kv.first) - 1, kv.second);

      // Allow defining a global node set containing all nodes.
      if (meshData.addGlobalNodeSet) {
        NEM::MSH::EXOMesh::ndeSetType ns;
        ns.id = ++ins;
        ns.nNde = em2.getNumberOfNodes();
        ns.name = meshData.addGlobalNodeSet.value();
        ns.ndeIdOffset = ndeIdOffset;
        for (int iNde = 0; iNde < ns.nNde; ++iNde)
          ns.ndeIds.emplace_back(iNde + 1);
        em2.addNdeSet(ns);
      }

      em->stitch(em2);
      ieb += em2.getNumberOfElementBlocks();
      ins += em2.getNumberOfNodeSets();
      iss += em2.getNumberOfSideSets();
      ndeIdOffset += em2.getNumberOfNodes();
      elmIdOffset += em2.getNumberOfElements();
      continue;
    }

    // reading input mesh
    meshBase *mb = meshBase::Create(meshData.meshFileName);
    // mb->write("exo_inp_mesh_"+std::to_string(iMsh)+".vtu");

    // adding information to exodusII object
    int ndeIdOffset_local = 0;
    int elmIdOffset_local = 0;

    ConversionDriver::genExo(
        mb, em, ndeIdOffset, elmIdOffset, ins, ieb, iss, meshData.meshName,
        meshData.usePhys, ndeIdOffset_local, elmIdOffset_local,
        meshData.makeFreeSurfSS, meshData.splitTopBotSS, meshData.sideSetNames);

    // offsetting starting ids for next file
    ndeIdOffset += ndeIdOffset_local;
    elmIdOffset += elmIdOffset_local;
    // clean up
    delete mb;
  }

  // writing the file
  em->write();
  em->report();

  // performing post-processing tasks
  if (this->opts_.needsPostProc) {
    for (int iTsk = 0; iTsk < this->opts_.numTasks; iTsk++) {
      const auto &ppFName = this->opts_.tasks[iTsk].taskFile;
      std::cout << "Reading Post Processing JSON file " << iTsk << std::endl;
      std::ifstream inputStream(ppFName);
      if (!inputStream.good() || nemAux::find_ext(ppFName) != ".json") {
        std::cerr << "Error opening file " << ppFName << std::endl;
        exit(1);
      }
      if (nemAux::find_ext(ppFName) != ".json") {
        std::cerr << "Input File must be in .json format" << std::endl;
        exit(1);
      }
      jsoncons::json ppJson;
      inputStream >> ppJson;
      ConversionDriver::procExo(ppJson, this->file_.outputFile, em);
    }

    // writing augmented exo file
    std::cout << "writing exodus file" << std::endl;
    em->write();
    em->report();
  }

  // clean up
  delete em;
}

}  // namespace DRV
}  // namespace NEM
