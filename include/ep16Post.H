#ifndef ep16Post_H
#define ep16Post_H

#ifdef HAVE_EPIC

// Nemosys headers
#include "nemosys_export.h"
#include <jsoncons/json.hpp>
#include "convexContainer.H"
#include "kmeans.H"
#include "convexContainer.H"
#include <memory>

// other
#include <iostream>

#define _kmeans_max_itr 100

namespace NEM {

namespace EPC {

/**
  @brief Implements methods for post-processing EPIC outputs.
**/
class NEMOSYS_EXPORT ep16Post
{
  public:
    /**
      @breif Contructs using input JSON object containing
      parameters
      @param inputjson is the input JSON objec
    **/
    explicit ep16Post(const jsoncons::json &inputjson):
        _jstrm(std::move(inputjson))
    {};

    // destructor
    ~ep16Post()
    {
      if (cont)
          delete cont;
      if (kmeans)
          delete kmeans;
    };

    /**
      @brief Static factory methods used to construct the object with given
      JSON object.
      @param jsonFName is the full path to the JSON input file
      @param inputjson is the JSON input object
      @param ret is the return value if successful it will be 0
    **/
    static ep16Post *readJSON(const std::string &jsonFName);
    static ep16Post *readJSON(const jsoncons::json &inputjson, int& ret);

  public:
    /**
      @breif Default JSON processing method.
      @return a non-zero value incase issues happen
    **/
    int readJSON();

    /**
      @brief Implements a combined KMeans-Quickhull method to identify convex
      objects from a point cloud
      @param jtsk is the JSON object describing the task parameters.
    **/
    int procErode(const jsoncons::json &jtsk);

  public:
    jsoncons::json _jstrm;
    NEM::MTH::KMeans* kmeans;
    NEM::GEO::convexContainer* cont;
};

} // namespace EPC 

} // namespace NEM

#endif // HAVE_EPIC

#endif // ep16Post_H
