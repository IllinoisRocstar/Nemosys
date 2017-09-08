#ifndef spheres_H
#define spheres_H

#include <iostream>
#include <fstream>
#include <ostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cstring>
#include <sstream>
#include <vector>
#include<string>

using std::cout; using std::endl;
using std::ofstream; using std::string;
using std::istream; using std::ostream;
using std::ifstream; using std::vector;
using std::size_t; using std::stringstream;

class sphere {
public:
  // constructors
  sphere(){}  
  sphere(double x_, double y_, double z_, double r_)
  : x(x_), y(y_), z(z_), r(r_) {}
  ~sphere(){}
  // access   
  double X() const { return x; }
  double Y() const { return y; }
  double Z() const { return z; }
  double R() const { return r; }

  // member
  bool in_sphere(std::vector<double> point);    

private:
  double x;
  double y;
  double z;
  double r;
};

typedef struct sphere_string sphere_string;
struct sphere_string 
{
  vector<sphere> spheres;
  vector<std::string> strings;
};

// read spheres from istream
sphere_string readSpheres(std::istream& inputStream);
// read spheres from file
sphere_string readSpheres(std::string filename);
// write spheres to outupt stream
void writeSpheres(const std::vector<sphere> s, ostream& outputStream);
#endif
