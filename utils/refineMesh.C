#include <RefineDriver.H>

int main(int argc, char* argv[])
{
  if (argc != 5)
  {
    exit(1);
  }
  
  std::string mesh(argv[1]);
  std::string method(argv[2]);
  std::string arrayName(argv[3]);
  std::string ofname(argv[4]);

  RefineDriver* driver = new RefineDriver(mesh, method, arrayName, 0.0, 1, 0.0,ofname); 
  if (driver) delete driver;
  return 0;

}
