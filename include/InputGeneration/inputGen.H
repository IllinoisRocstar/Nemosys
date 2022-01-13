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
#ifndef NEMOSYS_INPUTGEN_H_
#define NEMOSYS_INPUTGEN_H_

// headers
#include "nemosys_export.h"

#include <jsoncons/json.hpp>

#include <map>

class meshBase;

enum class inpFileType { INPGEN_TXT, INPGEN_CGNS, INPGEN_VTK };

// class declarations
// inputGen is an abstract class that provides a general-purpose,
// template structure for generating input needed for the
// high fidelity physics codes.
class NEMOSYS_EXPORT inputGen {
  // constructors and destructors
  public:
    explicit inputGen(jsoncons::json jsnStrm)
        : _jstrm(std::move(jsnStrm)), _key("generic") {};

    virtual ~inputGen() {
      for (auto &&inp : _inp)
        delete inp.second;
    }

  // basic methods
  public:
    // setting keyword for each input file
    void _set_key(const std::string &key) { _key = key; }

    virtual void setNameType(const std::string &fname, inpFileType ftyp,
                             const std::string &key = std::string());
    // order of writing lines into the file
    virtual void setOrder(const std::vector<std::string> &ord,
                          const std::string &key = std::string());
    virtual std::vector<std::string>
    getOrder(const std::string &key = std::string());
    virtual void pushOrder(const std::string &ord,
                           const std::string &key = std::string());
    // setting mesh database
    virtual void setMsh(meshBase *mb,
                        const std::string &key = std::string());
    // string used for commenting for each input file
    virtual void setCmntStr(const std::string &cmstr,
                            const std::string &key = std::string());
    virtual std::string getCmntStr(const std::string &key = std::string());
    // analyzing input and prepare for writing to the file
    virtual void process() = 0;
    // wrting into the file, implementation can be different
    virtual void write(const std::string &key = std::string()) const;
    // reading input file, implementation can be different
    virtual void read(const std::string &fname) = 0;
    // closing input file opened by read, implementation can be different
    virtual void close(const std::string &fname) = 0;

  // commonly used feature specifications
  public:
    // adding element sets and node sets
    virtual void addNdeSet() = 0;
    virtual void addElmSet() = 0;

  // attribute assignment
  public:
    virtual void addMat() = 0;
    virtual void addBC() = 0;
    virtual void addIC() = 0;
    virtual bool addMisc() = 0;

  // internal processing
  public:
    virtual void wrtPre(const std::string &_tsk, const std::string &__tsk) = 0;
    virtual void wrtCmnt(const std::string &cmnt) = 0;
    virtual void wrtMsh(const std::string &_tsk, const std::string &__tsk) = 0;
    virtual void wrtNdeCrds() = 0;
    virtual void wrtElmCon() = 0;
    virtual void wrtPost() = 0;
    virtual void wrtNdeSet() = 0;
    virtual void wrtElmSet() = 0;
    virtual void wrtMisc(const std::string &_tsk, const std::string &__tsk) = 0;
    virtual void edit(const std::string &_tsk, const std::string &__tsk) = 0;

  // raw methods
  protected:
    virtual void _write(std::stringstream &ss)
    {
      auto key_it = _inp.find(_key);
      if (key_it == _inp.end())
        _inp.insert(
            std::pair<std::string, std::stringstream *>(
                _key,
                new std::stringstream()
            )
        );

      *(_inp.at(_key)) << ss.str() << std::endl;
      ss.clear();
      ss.str(std::string());
    };

  protected:
    jsoncons::json _jstrm;
    std::string _key;
    std::map<std::string, std::stringstream *> _inp;
    std::map<std::string, std::vector<meshBase *>> _mb;
    std::map<std::string, std::vector<std::string>> _ord;
    std::map<std::string, std::string> _fn;
    std::map<std::string, inpFileType> _tpe;
    std::map<std::string, std::string> _cmnt;
};

#endif  // NEMOSYS_INPUTGEN_H_
