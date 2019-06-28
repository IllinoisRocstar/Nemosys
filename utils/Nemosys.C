#include "NemDriver.H"
#include "AuxiliaryFunctions.H"

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
    exit(1);
  }

  std::string fname = argv[1];
  std::ifstream inputStream(fname);
  if (!inputStream.good() || nemAux::find_ext(fname) != ".json")
  {
    std::cerr << "Error opening file " << fname << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  if (inputjson.is_array())
    for (const auto &prog : inputjson.array_range())
    {
      NemDriver *nemdrvobj = NemDriver::readJSON(prog);
      delete nemdrvobj;
    }
  else
  {
    NemDriver *nemdrvobj = NemDriver::readJSON(inputjson);
    delete nemdrvobj;
  }

  return 0;
}

//TODO
// also overload each driver constructor with one that takes pointer to resuling mesh
// with that, we can maybe do more things in memory

// add netgen uniform option to Refinement program (i.e. Refinement Method = "uniformMAd" || "uniformNG")
