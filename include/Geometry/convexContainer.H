/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_CONVEXCONTAINER_H_
#define NEMOSYS_CONVEXCONTAINER_H_

// Nemosys headers
#include "nemosys_export.h"

#include <iostream>
#include <iterator>
#include <fstream>
#include <cstddef>
#include <vector>
#include <map>
#include <cassert>
#include <cmath>

#include <QuickHull.hpp>
#include <Structs/Vector3.hpp>

// Vector struct used for convexContainer object
namespace NEM {

namespace MTH {

// TODO: move this to math elementaries

struct Vector {
  double x, y, z;

  Vector operator-(Vector p) const 
  {
    return Vector{x - p.x, y - p.y, z - p.z};
  }

  Vector cross(Vector p) const 
  {
    return Vector{
      y * p.z - p.y * z,
      z * p.x - p.z * x,
      x * p.y - p.x * y
    };
  }

  double dot(Vector p) const 
  {
    return x * p.x + y * p.y + z * p.z;
  }

  double norm() const 
  {
    return std::sqrt(x*x + y*y + z*z);
  }
};

} // namespace GEO

} // namespace NEM



namespace NEM {

namespace GEO {


struct Face {
  std::vector<NEM::MTH::Vector> v;

  NEM::MTH::Vector normal() const 
  {
    assert(v.size() > 2);
    NEM::MTH::Vector dir1 = v[1] - v[0];
    NEM::MTH::Vector dir2 = v[2] - v[0];
    NEM::MTH::Vector n  = dir1.cross(dir2);
    double d = n.norm();
    return NEM::MTH::Vector{n.x / d, n.y / d, n.z / d};
  }
};

/** @brief An implementation of convex container object for point
           The class implements methods to generate a convex hull for
           a point cloud, perform point location, and output to STL
           format.
    @note Only use for convex objects.
**/
class NEMOSYS_EXPORT convexContainer
{

  public:
  /** @brief default constructor 
  **/
  convexContainer();

  /** @brief construct from a vector point coordinates 
      @param verts vector of vector of doubles of point coordinats 
  **/
  convexContainer(std::vector<std::vector<double> >& verts);
  convexContainer(std::vector<quickhull::Vector3<double> >& verts);

  virtual ~convexContainer() {};

  /** @brief re-/setting vertice coordinates
      @param verts vector of vector of doubles of point coordinats 
  **/
  virtual void setVertex(std::vector<std::vector<double> >& verts);

  
  /** @brief calculating the convex hull
  **/
  virtual void computeConvexHull();

  
  /** @brief check for point containment
      @param p vector of coordinates of the point
  **/
  bool isInConvexPoly(NEM::MTH::Vector const& p);
  bool isInConvexPoly(const std::vector<double>& p);

  
  /** @brief generate STL triangulation from the convex hull
      @param file_name full path to the output STL file 
  **/
  void toSTL(std::string file_name) const;

  private:
  bool _isReady;
  std::vector<quickhull::Vector3<double> > vrts;
  quickhull::QuickHull<double>* qHull;
  std::vector<Face> fv;

};


} // namespace GEO

} // namespace NEM

#endif  // NEMOSYS_CONVEXCONTAINER_H_
