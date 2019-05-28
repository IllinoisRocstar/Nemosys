#ifndef STLTOVTK_H
#define STLTOVTK_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <functional>
#include <algorithm>

// Nemosys headers
#include <AuxiliaryFunctions.H>

class triangle
{
  public:
    triangle(std::vector<double> _p1, std::vector<double> _p2, std::vector<double> _p3):
      p1(_p1),p2(_p2),p3(_p3) {};
  
    ~triangle(){};   
  
  public:
    std::vector<double> p1;
    std::vector<double> p2;
    std::vector<double> p3;

};

void exportStlToVtk(std::string fname);

#endif
