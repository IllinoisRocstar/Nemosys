// -*- C++ -*-

#ifndef _H_MADVECTOR3
#define _H_MADVECTOR3

#ifdef _HAVE_GMSH_
  
#include "gmsh/SVector3.h"
namespace MAd {
  typedef SVector3 vec3;
}

#else
  
#include "MAdMessage.h"

namespace MAd {

  class MAdVector3
  {
  public:
    MAdVector3() {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Metric not implemented in MAdLib, Gmsh required");
    }
    MAdVector3(double x, double y, double z) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Metric not implemented in MAdLib, Gmsh required");
    }
    MAdVector3(double v) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Metric not implemented in MAdLib, Gmsh required");
    }
    MAdVector3(const double *array) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Metric not implemented in MAdLib, Gmsh required");
    }
    MAdVector3(const MAdVector3& v) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Metric not implemented in MAdLib, Gmsh required");
    }
    inline double x(void) const { return -1.; }
    inline double y(void) const { return -1.; }
    inline double z(void) const { return -1.; }
    inline double norm() { return -1.; }
    inline double normSq() { return -1.; }
    double normalize() { return -1.; }
    void negate() {}  
    double &operator[](int i){ }
    double operator[](int i) const { return -1.; }
    double &operator()(int i){ }
    double operator()(int i) const { return -1.; }
    MAdVector3 & operator += (const MAdVector3 &a) { return *this; }
    MAdVector3 & operator -= (const MAdVector3 &a) { return *this; }
    MAdVector3 & operator *= (const MAdVector3 &a) { return *this; }
    MAdVector3 & operator *= (const double v) { return *this; }
    MAdVector3 & operator = (double v) { return *this; }
    operator double *() { return NULL; }
    void print(std::string name="") const {}
  };
  
  typedef MAdVector3 vec3;
  
  inline MAdVector3 crossprod(const MAdVector3 &a, const MAdVector3 &b)
  { return MAdVector3(); }

}

#endif

#endif
