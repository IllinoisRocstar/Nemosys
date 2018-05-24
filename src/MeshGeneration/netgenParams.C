#include <netgenParams.H>


netgenParams::netgenParams()
{
  uselocalh                   = 1;
  maxh                        = 1000.0;
  fineness                    = 0.5;
  grading                     = 0.3;
  elementsperedge             = 2.0;
  elementspercurve            = 2.0; 
  closeedgeenable             = 0;
  closeedgefact               = 2.0;
  second_order                 = 0;
  meshsize_filename           = "null";
  quad_dominated              = 0;
  optvolmeshenable            = 1;
  optsteps_2d                 = 3;
  optsteps_3d                 = 3;
  invert_tets                 = 0;
  invert_trigs                = 0;
  check_overlap               = 1;
  check_overlapping_boundary  = 1;
  refine_with_geom            = 0;
  refine_without_geom         = 0;  
}

