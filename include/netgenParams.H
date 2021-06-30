#ifndef NETGENPARAMS_H
#define NETGENPARAMS_H

#include "nemosys_export.h"
#include "meshingParams.H"

#include <string>

class NEMOSYS_EXPORT netgenParams : public meshingParams
{
  public:
    // initializing params with defaults
    netgenParams() = default;

  public:
    bool uselocalh{true};
    double maxh{1000.0};
    double fineness{0.5};
    double grading{0.3};
    double elementsperedge{2.0};
    double elementspercurve{2.0};
    bool closeedgeenable{false};
    double closeedgefact{2.0};
    bool second_order{false};
    std::string meshsize_filename{};
    bool quad_dominated{false};
    bool optvolmeshenable{true};
    int optsteps_2d{3};
    int optsteps_3d{3};
    bool invert_tets{false};
    bool invert_trigs{false};
    bool check_overlap{true};
    bool check_overlapping_boundary{true};
    bool refine_with_geom{false};
    bool refine_without_geom{false};
};

#endif
