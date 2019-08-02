// Basic storage class for an n-dimensional point with a cluster assignment.
// Originally used for KMeans calculations
// Adopted from Felix Duvallet

#ifndef __KMEANS_POINT_H__
#define __KMEANS_POINT_H__

#include "nemosys_export.h"

#include <vector>
#include <iostream>

namespace NEM {

namespace GEO {

/** @brief Point class implements a simple point representation.
**/
class NEMOSYS_EXPORT Point {
  public:
  Point() {}

  /** @brief Initialize the number of dimensions, optionally set all values to zero
      @param num_dimensions number of dimensions
      @param init_zero if true initialize with zeros
  **/
  Point(int num_dimensions, bool init_zeros = true);

  /** @brief Initialize the point by taking the coordinates
      @param x coordinate in the x-direction
      @param y coordinate in the y-direction
      @param z coordinate in the z-direction
  **/
  Point(double x, double y, double z);

  /** @brief Initialize by vector
      @param vector coordinates
  **/
  Point(const std::vector<double> &vector);

  ~Point() { }

  // 

  /** @brief Compute distance between two points
      @param p1 first point
      @param p2 second point
  **/
  static double distance(const Point &p1, const Point &p2);

  /** @brief Adds a point to the current point
      @param point point object
  **/
  void add(const Point &point);

  /** @brief Update the cluster assignment. Returns true if cluster assignment
      changed, false if it stayed the same.
      @param k new cluster id
  **/
  bool update(int k);

  // Members: the data, the number of dimensions, and the cluster ID.
  std::vector<double> data_;
  int dimensions_;
  int cluster_;

  friend std::ostream &operator<<(std::ostream &target, const Point &point);
};

} // namespace GEO

} // namespace NEM

#endif  // __KMEANS_POINT_H_
