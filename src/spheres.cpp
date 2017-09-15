#include <iostream>
#include <fstream>
#include <ostream>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>
#include <spheres.H>

using std::cout; using std::endl;
using std::ofstream; using std::string;
using std::istream; using std::ostream;
using std::ifstream; using std::vector;
using std::size_t; using std::stringstream;

bool sphere::in_sphere(std::vector<double> point)
{
  if (point.size() != 3) {
    std::cerr << "Point must be triplet" << std::endl;
    exit(3);
  }
  double dist = pow(x-point[0],2) +
                pow(y-point[1],2) +
                pow(z-point[2],2);
  if (dist <= pow(r,2))
    return true;
  return false;
}


// read spheres from istream 
//automatically supports derived classes of istream by inheritance
sphere_string readSpheres(istream& inputStream)
{
  double tmp;
  string line;
  vector<sphere> spheres;
  vector<string> strings;
  sphere_string spherestring;
  int i = 1;
  // grab each line in file
  while (getline(inputStream, line)) {
    // only consider sphere lines
    if (line.find("Sphere(") != -1) {
      // defining bounds for actual data
      size_t beg = line.find("{");
      size_t end = line.find("}");
      if (beg != -1 && end != -1) {
        // creating substring from beg-end bounds of line
        string data = line.substr(beg+1, end-beg-1);
        // creating stream for substring to allow polymorphic output to var
        stringstream ss(data); 
        // declaring vector for sphere coordinates
        vector<double> coords;
        // pushing data from substring stream to tmp variable
        while (ss >> tmp) {
          coords.push_back(tmp);
          // if comma encountered, move on
          if (ss.peek() == ',')
            ss.ignore();
        }
        // populating sphere vector to be returned
        spheres.push_back(sphere(coords[0], coords[1], coords[2], coords[3]));      
      }
    }
    else if (line.find("Physical Volume") != -1) {
      size_t beg = line.find("(");
      size_t end = line.find(")");
      if (beg != -1 && end != -1) {
        string data = line.substr(beg+1, end-beg-1);
        size_t first = data.find_first_of("\"");  
        size_t last = data.find("\"", first+1);       
        size_t num1 = line.find("{");
        size_t num2 = line.find("}");
        if  (first != -1 && last != -1 && num1 != -1 && num2 != -1) { 
          string sph_num = line.substr(num1+1, num2-num1 -1);
          stringstream ss(sph_num);
          std::vector<double> inds;
          while (ss >> tmp) {
            inds.push_back(tmp);
            if (ss.peek() == ':')
              ss.ignore();
          }
          for (int i = 0; i < inds[1] - inds[0] +1; ++i)      
            strings.push_back(data.substr(first+1,last-first-1));
        }
        else {
          std::cout << "Error parsing .geo file at line " << i << std::endl;
          exit(1);
        } 
      }
      else {
        std::cout << "Error parsing .geo file at line " << i << std::endl;
        exit(1);
      }           
    }
    i+=1;
  }
  spherestring.spheres = spheres;
  spherestring.strings = strings;
  return spherestring;
}

// read spheres from file (using ifstreams inheritance of istream methods)
sphere_string readSpheres(string filename)
{
  ifstream inputStream(filename.c_str());
  if (!inputStream.good()) {
    cout << "Error: " << filename << " not found" << endl;
    exit(1);
  }
  return readSpheres(inputStream);
}

// write spheres to output stream
void writeSpheres(const std::vector<sphere> s, ostream& outputStream)
{
  for (int i = 0; i < s.size(); ++i) {
    outputStream << s[i].X() << " " << s[i].Y() << " " 
                 << s[i].Z() << " " << s[i].R() << std::endl; 
  }
}







