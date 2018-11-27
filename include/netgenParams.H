#ifndef NETGENPARAMS_H
#define NETGENPARAMS_H
#include <meshingParams.H>
#include <string>

class netgenParams : public meshingParams
{
  
  public:
    // initializing params with defaults
    netgenParams();
    ~netgenParams(){};
  
  public:
    bool uselocalh;
    double maxh;
    double fineness;  
    double grading;
    double elementsperedge;
    double elementspercurve;
    bool closeedgeenable;
    double closeedgefact;
    bool second_order;
    std::string meshsize_filename;
    bool quad_dominated;
    bool optvolmeshenable;
    int optsteps_2d;
    int optsteps_3d;
    bool invert_tets;
    bool invert_trigs;
    bool check_overlap;
    bool check_overlapping_boundary;
    bool refine_with_geom;
    bool refine_without_geom;
};



#endif
