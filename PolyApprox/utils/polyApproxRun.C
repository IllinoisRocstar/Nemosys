#include <polyApprox.H>
#include <fstream>
#include <sstream>
#include <iomanip>


int main(int argc, char* argv[])
{
  std::ifstream inputStream(argv[1]);
  std::ifstream inputStream1(argv[2]);
  std::string line,line1;
  VectorXd func(200);
  std::vector<std::vector<double>> coords(200);
  int i = 0;
  while (getline(inputStream, line) && getline(inputStream1, line1))
  {
    coords[i].resize(3);
    std::stringstream ss(line);
    ss >> func(i);
    std::stringstream ss1(line1);
    ss1 >> coords[i][0] >> coords[i][1] >> coords[i][2];
    i+=1;
  }
 

  for (int i = 0; i < 200; ++i)
  {
    std::cout << coords[i][0] << " " << coords[i][1] << " " << coords[i][2] << std::endl;
  }
 
  std::vector<double> node = {0.199088,0.0939398,0.582722};

  polyApprox polyapprox(1,coords);
  polyapprox.computeCoeff(func); 
  std::cout << polyapprox.eval(node) << std::endl; 

  return 0;
}
