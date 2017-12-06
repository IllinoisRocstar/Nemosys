#include <NemDrivers.H>


int main(int argc, char* argv[])
{

  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
    exit(1);
  }
  std::string fname = argv[1];
  RefineDriver* refdrvobj = RefineDriver::readJSON(fname);
  if (refdrvobj) 
  { 
    delete refdrvobj;
    refdrvobj = 0;
  }
  return 0;
}
