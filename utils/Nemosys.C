#include <NemDrivers.H>

int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
    exit(1);
  }
  std::string fname = argv[1];
  NemDriver* nemdrvobj = NemDriver::readJSON(fname);
  if (nemdrvobj) 
  { 
    delete nemdrvobj;
    nemdrvobj = 0;
  }
  return 0;

}
