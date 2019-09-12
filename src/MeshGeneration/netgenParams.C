#include "netgenParams.H"

netgenParams::netgenParams()
{
  uselocalh                   = true;
  maxh                        = 1000.0;
  fineness                    = 0.5;
  grading                     = 0.3;
  elementsperedge             = 2.0;
  elementspercurve            = 2.0;
  closeedgeenable             = false;
  closeedgefact               = 2.0;
  second_order                = false;
  meshsize_filename           = "null";
  quad_dominated              = false;
  optvolmeshenable            = true;
  optsteps_2d                 = 3;
  optsteps_3d                 = 3;
  invert_tets                 = false;
  invert_trigs                = false;
  check_overlap               = true;
  check_overlapping_boundary  = true;
  refine_with_geom            = false;
  refine_without_geom         = false;
}
