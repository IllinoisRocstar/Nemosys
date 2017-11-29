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

// flattens vector of vectors
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v);

// folds vector into vector of vectors
template <typename T>
std::vector<std::vector<T>> fold(const std::vector<T>& v, int dim);

// adds two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& x, 
                         const std::vector<T>& y);

// compute L2 norm of vec
inline double L2_Norm(const std::vector<double>& x);

// computes reciprocal of argument
template <typename T>
T reciprocal(T x);

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x);

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

// multiplies vector by scalar
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x);

// get average and stdev of values
inline std::vector<double> getMeanStdev(std::vector<double>& x);

// generates boolean array with 1 if value >= tol, 0 otherwise
inline std::vector<int> cellsToRefine(std::vector<double>& values, double tol);

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


// compute L2 norm of vec
double L2_Norm(const std::vector<double>& x)
{
  double result = 0.0;
  for (int i = 0; i < x.size(); ++i)
    result += x[i]*x[i];
  
  return std::sqrt(result);
}

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
  
  std::vector<T> result;
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), y.begin(), 
                 std::back_inserter(result), std::plus<T>());
  return result;
} 

// multiplies vector by scalar
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x)
{
  std::vector<T> result;
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), result.begin(), std::bind1st(std::multiplies<T>(),a));
  return result;
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
std::vector<double> getMeanStdev(std::vector<double>& x)
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
