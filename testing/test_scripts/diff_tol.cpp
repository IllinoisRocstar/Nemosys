#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>

int main (int argc, char* argv[])
{
  if (argc < 4) 
  {
    std::cout << "Usage: " << argv[0] << " file1 file2 TOL" << std::endl;
    exit(1);
  }

  std::ifstream inputstream1(argv[1]), inputstream2(argv[2]);

  if (!inputstream1.is_open() || !inputstream2.is_open()) 
  {
    std::cout << "Failed to open file" << std::endl;
    exit(2);
  }
  
  std::string line1, line2, tmp_str;
  double val1, val2, TOL;
  TOL = atof(argv[3]);
  int i = 0;
  while (getline(inputstream1, line1) && getline(inputstream2, line2)) 
  {
    ++i;
    std::istringstream in1(line1);
    std::istringstream in2(line2);
    while (in1 >> val1 && in2 >> val2) 
    {
      if (std::abs(val1-val2) > TOL) 
      {
        std::cout << "diff steps out of tolerance at line " << i << std::endl;
        exit(1);
      }
    }
  }
  return EXIT_SUCCESS;
}

