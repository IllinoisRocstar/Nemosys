#include <fstream>

#include "AuxiliaryFunctions.H"
#include "NemDriver.H"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
    exit(1);
  }

  std::string fname = argv[1];
  std::ifstream inputStream(fname);
  if (!inputStream.good() || nemAux::find_ext(fname) != ".json") {
    std::cerr << "Error opening file " << fname << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  if (inputjson.is_array())
    for (const auto &prog : inputjson.array_range()) {
      NEM::DRV::NemDriver *nemdrvobj = NEM::DRV::NemDriver::readJSON(prog);
      delete nemdrvobj;
    }
  else {
    NEM::DRV::NemDriver *nemdrvobj = NEM::DRV::NemDriver::readJSON(inputjson);
    delete nemdrvobj;
  }

  return 0;
}

// TODO: also overload each driver constructor with one that takes pointer to
//       resulting mesh with that, we can maybe do more things in memory.
//       Add Netgen uniform option to Refinement program (i.e. Refinement Method
//       = "uniformMAd" || "uniformNG")
