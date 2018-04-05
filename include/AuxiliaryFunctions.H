#ifndef AuxiliaryFunctions_H
#define AuxiliaryFunctions_H

#include <string>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <chrono>
//---------------------------Auxiliary Classes---------------------------------//

// class wrapping around chrono for timing methods
class Timer {
private:
  typedef std::chrono::time_point<std::chrono::system_clock> time_t;

public:
  Timer() : startTime(), stopTime() {}

  time_t start()   { return (startTime = std::chrono::system_clock::now()); }
  time_t stop()    { return (stopTime  = std::chrono::system_clock::now()); }
  double elapsed() { return std::chrono::duration_cast<std::chrono::milliseconds>
                                                      (stopTime-startTime).count(); }

private:
  time_t startTime, stopTime;
};

//----------------------------------------------------------------------//

//---------------------Auxiliarry Function Declarations---------------------------------//

/*
Reasons for inline functions:
 
  Per C++11 Standard, every program (after linking) shall contain exactly one definition 
  of every non-inline function or variable that is odr-used (one-definition-rule-used) in the
  program. An inline function shall be defined in every translation unit in which it 
  is odr-used.

  So, if functions are non-inline and defined in a header that is included in multiple 
  translation units (separate source files), the linker will see the same symbol is
  used more than once when merging object code. This violates the odr. 
  
  If functions are inline and defined in a header, an exception is made to the odr rule, 
  and it is defined in every translation unit in which it's used. Although inlining is a
  non-binding request, the inline keyword succeeds at telling the linker to tolerate
  multiple symbol definitions.
*/


// flattens vector of vectors
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v);

// folds vector into vector of vectors
template <typename T>
std::vector<std::vector<T>> fold(const std::vector<T>& v, int dim);

namespace nemAux {
// adds two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& x, 
                         const std::vector<T>& y);
// sub two vectors
template <typename T>
std::vector<T> operator-(const std::vector<T>& x,
                         const std::vector<T>& y);
// multiplies vector by scalar
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x);

// hadamard product of vectors
template <typename T>
std::vector<T> hadamard(const std::vector<T>& x, 
                        const std::vector<T>& y);
}

// compute 2 norm of vec
inline double l2_Norm(const std::vector<double>& x);

// computes reciprocal of argument
template <typename T>
T reciprocal(T x);

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x);

// checks if vector has zeros
inline bool hasZero(std::vector<double>& x);


// find minmax of vector, excluding inf 
inline std::vector<double> getMinMax(const std::vector<double>& x);

// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
inline double scale_to_range(double x, const std::vector<double>& xminmax, 
                      const std::vector<double>& yminmax);

// scales each number in vector in range [xmin, xmax] to [ymin,ymax]
inline void scale_vec_to_range(std::vector<double>& x, 
                        const std::vector<double>& xminmax,
                        const std::vector<double>& yminmax);

// get average and stdev of values
inline std::vector<double> getMeanStdev(const std::vector<double>& x);

// generates boolean array with 1 if value >= tol, 0 otherwise
inline std::vector<int> cellsToRefine(std::vector<double>& values, double tol);

// generates boolean array with 1 if value >= mean+/-dev, 0 otherwise
inline std::vector<int> cellsToRefineStdev(std::vector<double>& values, double mean, double dev);

// return selects cells if their value is above some threshold of the maximum
inline std::vector<int> cellsToRefineMaxdev(std::vector<double>& values, double dev);

// string trimming for consistent file names 
inline std::string trim_fname(std::string name, std::string ext);
// find file extension
inline std::string find_ext(std::string fname);

// print a vector
template<typename T>
void printVec(const std::vector<T>& v);
//-------------------------------------------------------------------------//


// --------------------------- Auxilliary Function Implementations -------------------------//

// flattens vector of vectors
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v)
{
  std::size_t size = 0;
  for (const auto& sub : v)
    size+= sub.size();
  std::vector<T> result;
  result.reserve(size);
  for (const auto& sub : v)
    result.insert(result.end(), sub.begin(), sub.end());
  return result;
}

// folds vector into vector of vectors, each of length dim
template <typename T>
std::vector<std::vector<T>> fold(const std::vector<T>& v, int dim)
{
  std::size_t size = v.size();
  if (size % dim != 0)
  {
    std::cout << "size must be divisible by dim" << std::endl;
    exit(1);
  }
  std::vector<std::vector<T>> result;
  result.resize(size/dim);
  std::cout << result.size() << std::endl;
  for (int i = 0; i < size/dim; ++i)
  {
    result[i].resize(dim);
    for (int j = 0; j < dim; ++j)
      result[i][j] = v[i*dim + j];
  }
  return result;
}


// compute 2 norm of vec
double l2_Norm(const std::vector<double>& x)
{
  double result = 0.0;
  for (int i = 0; i < x.size(); ++i)
    result += x[i]*x[i];
  
  return std::sqrt(result);
}

namespace nemAux {
// adds two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& x, 
                         const std::vector<T>& y)
{
  if (x.size() != y.size())
  {
    std::cout << "vectors must be same length for addition" << std::endl;
    exit(1);
  }
  
  std::vector<T> result(x.size());
  for (int i = 0; i < x.size(); ++i)
    result[i] = x[i] + y[i];
  /*result.reserve(x.size());
  std::transform(x.begin(), x.end(), y.begin(), 
                 std::back_inserter(result), std::plus<T>());
  */
  return result;
} 

// sub two vectors
template <typename T>
std::vector<T> operator-(const std::vector<T>& x,
                         const std::vector<T>& y)
{
  if (x.size() != y.size())
  {
    std::cout << "vectors must be same length for subtraction" << std::endl;
    exit(1);
  }
  
  std::vector<T> result(x.size());
  for (int i = 0; i < x.size(); ++i)
    result[i] = x[i] - y[i];
  /*result.reserve(x.size());
  std::transform(x.begin(),x.end(),y.begin(),
                 std::back_inserter(result), std::minus<T>());
  */
  return result;
}

// multiplies vector by scalar
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x)
{
  std::vector<T> result(x.size());
  for (int i = 0; i < x.size(); ++i)
    result[i] = a*x[i];
  //result.reserve(x.size());
  //std::transform(x.begin(), x.end(), result.begin(), std::bind1st(std::multiplies<T>(),a));
  return result;
}

template <typename T>
std::vector<T> hadamard(const std::vector<T>& x, 
                        const std::vector<T>& y)
{
  if (x.size() != y.size())
  {
    std::cout << "vectors must be same length for subtraction" << std::endl;
    exit(1);
  }
  
  std::vector<T> result(x.size());
  for (int i = 0; i < x.size(); ++i)
    result[i] = x[i]*y[i];
  return result; 
}
}

// computes reciprocal of argument
template <typename T>
T reciprocal(T x)
{
  return (T) 1/x ;
}

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x)
{
  std::transform(x.begin(), x.end(), x.begin(), reciprocal<T>);
}

bool hasZero(std::vector<double>& x)
{
    for (auto it=x.begin(); it!=x.end(); it++)
        if (*it==0.)
            return true;
    return false;
}

// find minmax of vector, excluding inf 
std::vector<double> getMinMax(const std::vector<double>& x)
{
  std::vector<double> result(2);
  std::vector<double> tmp;
  for (int i = 0; i < x.size(); ++i)
  {
    if (!std::isinf(x[i])) // exclude inf
    {
      tmp.push_back(x[i]); 
    }
  }
  auto minmax = std::minmax_element(tmp.begin(), tmp.end());
  result[0] = *minmax.first;
  result[1] = *minmax.second;
  return result;
} 

// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
double scale_to_range(double x, const std::vector<double>& xminmax, 
                      const std::vector<double>& yminmax)
{
  if (std::isinf(x))
    return yminmax[1];
  return yminmax[0] + (yminmax[1] - yminmax[0])*(x - xminmax[0])/(xminmax[1] - xminmax[0]);
}

// scales each number in vector in range [xmin, xmax] to [ymin,ymax]
void scale_vec_to_range(std::vector<double>& x, 
                        const std::vector<double>& xminmax,
                        const std::vector<double>& yminmax)
{
  for (int i = 0; i < x.size(); ++i)
    x[i] = scale_to_range(x[i], xminmax, yminmax);
}


// get average and stdev of values
std::vector<double> getMeanStdev(const std::vector<double>& x)
{
  std::vector<double> result(2);
  double ave = 0;
  for (int i = 0; i < x.size(); ++i)
    ave += x[i];
  ave /= x.size();

  double stdev = 0;
  for (int i = 0; i < x.size(); ++i)
    stdev += (x[i] - ave)*(x[i] - ave);
  stdev = std::sqrt(stdev/x.size());
  result[0] = ave;
  result[1] = stdev;
  return result;
}

// generates boolean array with 1 if value >= tol, 0 otherwise
std::vector<int> cellsToRefine(std::vector<double>& values, double tol)
{
  std::vector<int> result(values.size(),0);  
  for (int i = 0; i < values.size(); ++i)
  {
    if (values[i] > tol)
      result[i] = 1;
  }
  return result;
}

// generates boolean array with 1 if value >= mean+/-dev, 0 otherwise
std::vector<int> cellsToRefineStdev(std::vector<double>& values, double mean, double dev)
{
  std::cout << __FILE__ << __LINE__ << std::endl;
  if (std::abs(dev/mean) >= 1.0)
  {
    dev = 0.99*mean;
  }
  double th = mean + dev;
  double tl = mean - dev;
  std::vector<int> result(values.size(),0);  
  for (int i = 0; i < values.size(); ++i)
  {
    if (values[i] > th || values[i] < tl)
      result[i] = 1;
  }
  return result;
}

// generates boolean array with 1 if value >= mean+/-dev, 0 otherwise
std::vector<int> cellsToRefineMaxdev(std::vector<double>& values, double dev)
{
  dev = std::abs(dev);
  if (dev > 1.)
      dev = 1.;
  double max = *std::max_element(values.begin(), values.end());
  double hl = (1-dev)*max;
  std::vector<int> result(values.size(),0);  
  for (int i = 0; i < values.size(); ++i)
  {
    if (values[i] > hl)
      result[i] = 1;
  }
  return result;
}

// string trimming for consistent file names 
std::string trim_fname(std::string fname, std::string ext)
{
  size_t beg = 0;
  size_t end = fname.find_last_of('.');
  std::string name;
  if (end != -1)
  {
    name = fname.substr(beg,end);
    name.append(ext);
    return name;
  }

  else 
  {
    std::cout << "Error finding file extension for " << fname << std::endl;
    exit(1);
  }
}  

// find file extension
inline std::string find_ext(std::string fname)
{
  size_t last = fname.find_last_of('.');
  if (last != -1)
  {
    return fname.substr(last);  
  }
  else
  {
    std::cout << "Error finding file extension for " << fname << std::endl;
    exit(1);
  }

}
// print vector
template<typename T>
void printVec(const std::vector<T>& v)
{
  for (int i = 0; i < v.size(); ++i)
    std::cout << v[i] << " ";
  std::cout << std::endl;
}
#endif
