#ifndef STLTOVTK_H
#define STLTOVTK_H

#include "nemosys_export.h"

#include <string>
#include <vector>

class triangle
{
  public:
    triangle(std::vector<double> _p1,
             std::vector<double> _p2,
             std::vector<double> _p3) :
        p1(_p1), p2(_p2), p3(_p3) {};

    ~triangle() {};

  public:
    std::vector<double> p1;
    std::vector<double> p2;
    std::vector<double> p3;

};

void exportStlToVtk(std::string fname);

#endif
