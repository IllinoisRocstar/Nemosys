#ifndef spheres_H
#define spheres_H

#include <vector>
#include <string>
#include <istream>
#include <ostream>

#include "interp_export.h"

class INTERP_EXPORT sphere {
  public:
    // constructors
    sphere() = default;
    sphere(double x_, double y_, double z_, double r_)
        : x(x_), y(y_), z(z_), r(r_) {}
    ~sphere() = default;
    // access
    double X() const { return x; }
    double Y() const { return y; }
    double Z() const { return z; }
    double R() const { return r; }

    // member
    bool in_sphere(const std::vector<double> &point) const;

  private:
    double x;
    double y;
    double z;
    double r;
};

typedef struct sphere_string sphere_string;
struct INTERP_EXPORT sphere_string {
  std::vector<sphere> spheres;
  std::vector<std::string> strings;
};

// read spheres from istream
INTERP_EXPORT sphere_string readSpheres(const std::istream &inputStream);
// read spheres from file
INTERP_EXPORT sphere_string readSpheres(const std::string &filename);
// write spheres to output stream
INTERP_EXPORT void writeSpheres(const std::vector<sphere> &s,
                                 std::ostream &outputStream);

#endif // spheres_H
