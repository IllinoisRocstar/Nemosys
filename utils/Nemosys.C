#include <NemDrivers.H>

int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
    exit(1);
  }
  std::string fname = argv[1];
  std::ifstream inputStream(fname);
  if (!inputStream.good() || find_ext(fname) != ".json")
  {
    std::cout << "Error opening file " << fname << std::endl;
    exit(1);
  }
  if (find_ext(fname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  json inputjson;
  inputStream >> inputjson;
	for(const auto& prog : inputjson.array_range())
	{
			//TODO
			// also overload each driver constructor with one that takes pointer to resuling mesh
			// with that, we can maybe do more things in memory

			// add netgen uniform option to Refinement program (i.e. Refinement Method = "uniformMAd" || "uniformNG")
  		NemDriver* nemdrvobj = NemDriver::readJSON(prog);
  		if (nemdrvobj) 
  		{ 
  		  delete nemdrvobj;
  		  nemdrvobj = 0;
  		}
	}
  return 0;

}
