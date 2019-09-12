#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include <spheres.H>

bool sphere::in_sphere(const std::vector<double> &point) const
{
  if (point.size() != 3) {
    std::cerr << "Point must be triplet" << std::endl;
    exit(3);
  }
  double dist = pow(x - point[0], 2) +
                pow(y - point[1], 2) +
                pow(z - point[2], 2);
  return dist <= pow(r, 2);
}

// read spheres from istream
// automatically supports derived classes of istream by inheritance
sphere_string readSpheres(std::istream &inputStream)
{
  double tmp;
  std::string line;
  std::vector<sphere> spheres;
  std::vector<std::string> strings;
  sphere_string sphereString;
  int i = 1;
  // grab each line in file
  while (getline(inputStream, line)) {
    // only consider sphere lines
    if (line.find("Sphere(") != std::string::npos) {
      // defining bounds for actual data
      size_t beg = line.find('{');
      size_t end = line.find('}');
      if (beg != std::string::npos && end != std::string::npos) {
        // creating substring from beg-end bounds of line
        std::string data = line.substr(beg + 1, end - beg - 1);
        // creating stream for substring to allow polymorphic output to var
        std::stringstream ss(data);
        // declaring vector for sphere coordinates
        std::vector<double> coords;
        // pushing data from substring stream to tmp variable
        while (ss >> tmp) {
          coords.push_back(tmp);
          // if comma encountered, move on
          if (ss.peek() == ',')
            ss.ignore();
        }
        // populating sphere vector to be returned
        spheres.emplace_back(coords[0], coords[1], coords[2], coords[3]);
      }
    } else if (line.find("Physical Volume") != std::string::npos) {
      size_t beg = line.find('(');
      size_t end = line.find(')');
      if (beg != std::string::npos && end != std::string::npos) {
        std::string data = line.substr(beg + 1, end - beg - 1);
        size_t first = data.find_first_of('"');
        size_t last = data.find('"', first + 1);
        size_t num1 = line.find('{');
        size_t num2 = line.find('}');
        if (first != std::string::npos && last != std::string::npos
            && num1 != std::string::npos && num2 != std::string::npos) {
          std::string sph_num = line.substr(num1 + 1, num2 - num1 - 1);
          std::stringstream ss(sph_num);
          std::vector<double> inds;
          while (ss >> tmp) {
            inds.push_back(tmp);
            if (ss.peek() == ':')
              ss.ignore();
          }
          for (int i = 0; i < inds[1] - inds[0] + 1; ++i)
            strings.push_back(data.substr(first + 1, last - first - 1));
        } else {
          std::cout << "Error parsing .geo file at line " << i << std::endl;
          exit(1);
        }
      } else {
        std::cout << "Error parsing .geo file at line " << i << std::endl;
        exit(1);
      }
    }
    i += 1;
  }
  sphereString.spheres = spheres;
  sphereString.strings = strings;
  return sphereString;
}

// read spheres from file (using ifstreams inheritance of istream methods)
sphere_string readSpheres(const std::string &filename)
{
  std::ifstream inputStream(filename);
  if (!inputStream.good()) {
    std::cout << "Error: " << filename << " not found" << std::endl;
    exit(1);
  }
  return readSpheres(inputStream);
}

// write spheres to output stream
void writeSpheres(const std::vector<sphere> &s,
                  std::ostream &outputStream)
{
  for (const auto &i : s) {
    outputStream << i.X() << " " << i.Y() << " "
                 << i.Z() << " " << i.R() << std::endl;
  }
}
