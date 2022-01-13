/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_EP16POST_H_
#define NEMOSYS_EP16POST_H_

#ifdef HAVE_EPIC

#include "nemosys_export.h"

#include <memory>
#include <iostream>

#include <jsoncons/json.hpp>

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

  private:
    static constexpr int _kmeans_max_itr = 100;
    jsoncons::json _jstrm;
};

} // namespace EPC 

} // namespace NEM

#endif // HAVE_EPIC

#endif  // NEMOSYS_EP16POST_H_
