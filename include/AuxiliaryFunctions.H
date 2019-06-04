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
#include <iomanip>
#include <string.h>
#include <glob.h>

//--------------------------Start nemAux namespace----------------------------//
namespace nemAux {

//----------------------------Auxiliary Classes-------------------------------//

// class wrapping around std::chrono for timing methods
class Timer {
  private:
    typedef std::chrono::time_point<std::chrono::system_clock> time_t;

  public:
    Timer() : startTime(), stopTime() {}

    time_t start() { return startTime = std::chrono::system_clock::now(); }

    time_t stop() { return stopTime = std::chrono::system_clock::now(); }

    double elapsed()
    {
      return std::chrono::duration_cast<std::chrono::milliseconds>
          (stopTime - startTime).count();
    }

  private:
    time_t startTime, stopTime;
};

//----------------------------------------------------------------------------//

//----------------------Auxiliary Function Declarations-----------------------//

/*
Reasons for inline functions:
 
  Per C++11 Standard, every program (after linking) shall contain exactly one
  definition of every non-inline function or variable that is odr-used
  (one-definition-rule-used) in the program. An inline function shall be defined
  in every translation unit in which it is odr-used.

  So, if functions are non-inline and defined in a header that is included in
  multiple translation units (separate source files), the linker will see the
  same symbol is used more than once when merging object code. This violates the
  odr.
  
  If functions are inline and defined in a header, an exception is made to the
  odr rule, and it is defined in every translation unit in which it's used.
  Although inlining is a non-binding request, the inline keyword succeeds at
  telling the linker to tolerate multiple symbol definitions.
*/

// flattens vector of vectors
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &v);

// folds vector into vector of vectors
template<typename T>
std::vector<std::vector<T>> fold(const std::vector<T> &v, int dim);

// adds two vectors
template<typename T>
std::vector<T> operator+(const std::vector<T> &x,
                         const std::vector<T> &y);
// sub two vectors
template<typename T>
std::vector<T> operator-(const std::vector<T> &x,
                         const std::vector<T> &y);

// multiplies vector by scalar
template<typename T>
std::vector<T> operator*(T a,
                         const std::vector<T> &x);

// Hadamard product of vectors
template<typename T>
std::vector<T> hadamard(const std::vector<T> &x,
                        const std::vector<T> &y);

// search for pattern and return vector of matches (wrapper around posix glob)
inline std::vector<std::string> glob(const std::string &pattern);

// compute 2 norm of vec
template<typename T>
inline T l2_Norm(const std::vector<T> &x);
// computes reciprocal of argument
template<typename T>
T reciprocal(T x);
// IN PLACE, computes reciprocal of vector elements
template<typename T>
void reciprocal_vec(std::vector<T> &x);
// checks if vector has zeros
template<typename T>
inline bool hasZero(const std::vector<T> &x);
// find minmax of vector, excluding inf
template<typename T>
inline std::vector<T> getMinMax(const std::vector<T> &x);
// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
template<typename T>
inline T scale_to_range(T x,
                        const std::vector<T> &xminmax,
                        const std::vector<T> &yminmax);
// IN PLACE, scales each number in vector in range [xmin, xmax] to [ymin, ymax]
template<typename T>
inline void scale_vec_to_range(std::vector<T> &x,
                               const std::vector<T> &xminmax,
                               const std::vector<T> &yminmax);
// get average and stdev of values
template<typename T>
inline std::vector<T> getMeanStdev(const std::vector<T> &x);
// generates boolean array with 1 if value >= tol, 0 otherwise
template<typename T>
inline std::vector<bool> cellsToRefine(const std::vector<T> &values, T tol);
// generates boolean array with 1 if value >= mean+/-dev, 0 otherwise
template<typename T>
inline std::vector<bool>
cellsToRefineStdev(const std::vector<T> &values, T mean, T dev);
// return selects cells if their value is above some threshold of the maximum
template<typename T>
inline std::vector<bool>
cellsToRefineMaxdev(const std::vector<T> &values, T dev);
// string trimming for consistent file names
inline std::string trim_fname(const std::string &name, const std::string &ext);
// find file extension
inline std::string find_ext(const std::string &fname);
// IN PLACE, to lower
inline void toLower(std::string &str);
// IN PLACE, to upper
inline void toUpper(std::string &str);
// return from beginning to the location of str
inline std::string findToStr(const std::string &str, const std::string &ptrn);
// return from the location of str
inline std::string findFromStr(const std::string &str, const std::string &ptrn);
// returns current time as string
inline std::string getTimeStr();
// find file extension
inline std::string find_name(const std::string &fname);
// print a vector
template<typename T>
void printVec(const std::vector<T> &v);
// find if in vector
template<typename T>
bool valInVec(const std::vector<T> &v, T val);
// flip a pair
template<typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B> &p);
// flip a map into a multimap
template<typename A, typename B>
std::multimap<B, A> flip_map(const std::map<A, B> &src);
// check if coordinate is inside a bounding box
template<typename T>
bool isInBBox(const std::vector<T> &crd, const std::vector<T> &bb);
// get a vector of the keys from a map (which are sorted)
template<typename A, typename B>
std::vector<A> getSortedKeys(const std::map<A, B> &mapObj);
//----------------------------------------------------------------------------//


//-------------------Auxiliary Function Implementations-----------------------//

// flattens vector of vectors
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &v)
{
  std::size_t size = 0;
  for (const auto &sub : v)
    size += sub.size();

  std::vector<T> result;
  result.reserve(size);
  for (const auto &sub : v)
    result.insert(result.end(), sub.begin(), sub.end());
  return result;
}

// folds vector into vector of vectors, each of length dim
template<typename T>
std::vector<std::vector<T>> fold(const std::vector<T> &v, int dim)
{
  std::size_t size = v.size();
  if (size % dim != 0) {
    std::cout << "Size must be divisible by dim for folding" << std::endl;
    exit(1);
  }

  std::vector<std::vector<T>> result;
  result.resize(size / dim);
  for (std::size_t i = 0; i < size / dim; ++i) {
    result[i].resize(dim);
    for (std::size_t j = 0; j < dim; ++j)
      result[i][j] = v[i * dim + j];
  }
  return result;
}

// compute 2 norm of vec
template<typename T>
T l2_Norm(const std::vector<T> &x)
{
  T result = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i)
    result += x[i] * x[i];

  return std::sqrt(result);
}

// adds two vectors
template<typename T>
std::vector<T> operator+(const std::vector<T> &x,
                         const std::vector<T> &y)
{
  if (x.size() != y.size()) {
    std::cout << "Vectors must be same length for addition" << std::endl;
    exit(1);
  }

  std::vector<T> result(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    result[i] = x[i] + y[i];
/*
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), y.begin(),
                 std::back_inserter(result), std::plus<T>());
*/
  return result;
}

// sub two vectors
template<typename T>
std::vector<T> operator-(const std::vector<T> &x,
                         const std::vector<T> &y)
{
  if (x.size() != y.size()) {
    std::cout << "Vectors must be same length for subtraction" << std::endl;
    exit(1);
  }

  std::vector<T> result(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    result[i] = x[i] - y[i];
/*
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), y.begin(),
                 std::back_inserter(result), std::minus<T>());
*/
  return result;
}

// multiplies vector by scalar
template<typename T>
std::vector<T> operator*(T a,
                         const std::vector<T> &x)
{
  std::vector<T> result(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    result[i] = a * x[i];
/*
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), result.begin(),
                 std::bind1st(std::multiplies<T>(), a));
*/
  return result;
}

template<typename T>
std::vector<T> hadamard(const std::vector<T> &x,
                        const std::vector<T> &y)
{
  if (x.size() != y.size()) {
    std::cout << "Vectors must be same length for Hadamard product" << std::endl;
    exit(1);
  }

  std::vector<T> result(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    result[i] = x[i] * y[i];
  return result;
}

// computes reciprocal of argument
template<typename T>
T reciprocal(T x)
{
  return (T) 1 / x;
}

// IN PLACE, computes reciprocal of vector elements
template<typename T>
void reciprocal_vec(std::vector<T> &x)
{
  std::transform(x.begin(), x.end(), x.begin(), reciprocal<T>);
}

template<typename T>
bool hasZero(const std::vector<T> &x)
{
  for (auto it = x.begin(); it != x.end(); it++)
    if (*it == (T) 0)
      return true;
  return false;
}

// find minmax of vector, excluding inf
template<typename T>
std::vector<T> getMinMax(const std::vector<T> &x)
{
  std::vector<T> tmp;
  for (std::size_t i = 0; i < x.size(); ++i)
    if (!std::isinf(x[i])) // exclude inf
      tmp.push_back(x[i]);

  auto minmax = std::minmax_element(tmp.begin(), tmp.end());

  std::vector<T> result(2);
  result[0] = *minmax.first;
  result[1] = *minmax.second;
  return result;
}

// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
template<typename T>
T scale_to_range(T x,
                 const std::vector<T> &xminmax,
                 const std::vector<T> &yminmax)
{
  if (std::isinf(x))
    return yminmax[1];
  return yminmax[0] + (yminmax[1] - yminmax[0]) * (x - xminmax[0]) /
                      (xminmax[1] - xminmax[0]);
}

// IN PLACE, scales each number in vector in range [xmin, xmax] to [ymin,ymax]
template<typename T>
void scale_vec_to_range(std::vector<T> &x,
                        const std::vector<T> &xminmax,
                        const std::vector<T> &yminmax)
{
  for (std::size_t i = 0; i < x.size(); ++i)
    x[i] = scale_to_range(x[i], xminmax, yminmax);
}

// get average and stdev of values
template<typename T>
std::vector<T> getMeanStdev(const std::vector<T> &x)
{
  T ave = 0;
  for (std::size_t i = 0; i < x.size(); ++i)
    ave += x[i];
  ave /= x.size();

  T stdev = 0;
  for (std::size_t i = 0; i < x.size(); ++i)
    stdev += (x[i] - ave) * (x[i] - ave);
  stdev = std::sqrt(stdev / x.size());

  std::vector<T> result(2);
  result[0] = ave;
  result[1] = stdev;
  return result;
}

// generates boolean array with 1 if value >= tol, 0 otherwise
template<typename T>
std::vector<bool> cellsToRefine(const std::vector<T> &values, T tol)
{
  std::vector<bool> result(values.size(), false);
  for (std::size_t i = 0; i < values.size(); ++i)
    if (values[i] > tol)
      result[i] = true;
  return result;
}

// generates boolean array with 1 if value >= mean+/-dev, 0 otherwise
template<typename T>
std::vector<bool>
cellsToRefineStdev(const std::vector<T> &values, T mean, T dev)
{
  if (std::abs(dev / mean) >= (T) 1.0)
    dev = (T) 0.99 * mean;

  T th = mean + dev;
  T tl = mean - dev;

  std::vector<bool> result(values.size(), false);
  for (std::size_t i = 0; i < values.size(); ++i)
    if (values[i] > th || values[i] < tl)
      result[i] = true;
  return result;
}

// generates boolean array with 1 if value >= mean+/-dev, 0 otherwise
template<typename T>
std::vector<bool> cellsToRefineMaxdev(const std::vector<T> &values, T dev)
{
  dev = std::abs(dev);
  if (dev > (T) 1)
    dev = (T) 1;

  T max = *std::max_element(values.begin(), values.end());
  T hl = (1 - dev) * max;

  std::vector<bool> result(values.size(), false);
  for (std::size_t i = 0; i < values.size(); ++i)
    if (values[i] > hl)
      result[i] = true;
  return result;
}

// string trimming for consistent file names 
std::string trim_fname(const std::string &fname, const std::string &ext)
{
  size_t end = fname.find_last_of('.');
  if (end != std::string::npos) {
    std::string name = fname.substr(0, end);
    return name.append(ext);
  } else {
    std::cout << "Error finding file extension for " << fname << std::endl;
    exit(1);
  }
}

// find file extension
inline std::string find_ext(const std::string &fname)
{
  size_t last = fname.find_last_of('.');
  if (last != std::string::npos) {
    return fname.substr(last);
  } else {
    std::cout << "Error finding file extension for " << fname << std::endl;
    exit(1);
  }
}

// IN PLACE, to lower
inline void toLower(std::string &str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

// IN PLACE, to upper
inline void toUpper(std::string &str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

// return from beginning to the location of str
inline std::string findToStr(const std::string &str, const std::string &ptrn)
{
  size_t loc = str.find_first_of(ptrn);
  return str.substr(0, loc);
}

// return from the location of str
inline std::string findFromStr(const std::string &str, const std::string &ptrn)
{
  size_t loc = str.find_first_of(ptrn);
  return str.substr(loc + 1, std::string::npos);
}

// returns current time as string
inline std::string getTimeStr()
{
  time_t rawTime;
  struct tm *timeInfo;
  char buffer[80];

  time(&rawTime);
  timeInfo = localtime(&rawTime);

  strftime(buffer, sizeof(buffer), "%d-%m-%Y %I:%M:%S %p", timeInfo);
  std::string str(buffer);

  return str;
}

// find file name
inline std::string find_name(const std::string &fname)
{
  size_t first = fname.find_last_of('/');
  size_t last = fname.find_last_of('.');
  if (first != std::string::npos && last != std::string::npos) {
    return fname.substr(first + 1, last);
  } else {
    std::cout << "error finding file extension for " << fname << std::endl;
    exit(1);
  }
}

// print vector
template<typename T>
void printVec(const std::vector<T> &v)
{
  for (std::size_t i = 0; i < v.size(); ++i)
    std::cout << std::setprecision(15) << v[i] << " ";
  std::cout << std::endl;
}

// find is val is in vector
template<typename T>
bool valInVec(const std::vector<T> &v, T val)
{
  auto it = v.begin();
  for (; it != v.end(); ++it)
    if (*it == val)
      break;
  return it != v.end();
}

template<typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B> &p)
{
  return std::pair<B, A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B, A> flip_map(const std::map<A, B> &src)
{
  std::multimap<B, A> dst;
  std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                 flip_pair<A, B>);
  return dst;
}

template<typename T>
bool isInBBox(const std::vector<T> &crd, const std::vector<T> &bb)
{
  bool isIn = true;
  for (std::size_t i = 0; i < crd.size(); ++i) {
    if (!isIn)
      continue;
    if (crd[i] < bb[2 * i] || crd[i] > bb[2 * i + 1])
      isIn = false;
  }
  return isIn;
}

template<typename A, typename B>
std::vector<A> getSortedKeys(const std::map<A, B> &mapObj)
{
  if (!mapObj.empty()) {
    auto it = mapObj.begin();
    std::vector<A> sortedKeys(mapObj.size());
    int i = 0;
    while (it != mapObj.end()) {
      sortedKeys[i] = it->first;
      ++it;
      ++i;
    }
    return sortedKeys;
  } else {
    std::cerr << "Map is empty! No sorted keys to return." << std::endl;
    exit(1);
  }
}

// search for pattern and return vector of matches
std::vector<std::string> glob(const std::string &pattern)
{
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));
  int ret = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  if (ret == GLOB_NOSPACE || ret == GLOB_ABORTED) {
    globfree(&glob_result);
    std::cerr << "glob() failed with return value " << ret << std::endl;
    exit(1);
  } else if (ret == GLOB_NOMATCH) {
    return std::vector<std::string>();
  }

  std::vector<std::string> fnames;
  for (std::size_t i = 0; i < glob_result.gl_pathc; ++i) {
    fnames.emplace_back(glob_result.gl_pathv[i]);
  }

  globfree(&glob_result);
  return fnames;
}

//---------------------------End nemAux namespace-----------------------------//
} // namespace nemAux

#endif  // AuxiliaryFunctions_H
