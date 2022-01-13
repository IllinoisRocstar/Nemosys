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
#ifndef NEMOSYS_EP16PREP_H_
#define NEMOSYS_EP16PREP_H_

#ifdef HAVE_EPIC

// Nemosys headers
#include "nemosys_export.h"
#include "InputGeneration/inputGen.H"
#include "Mesh/exoMesh.H"

// other
#include <iostream>


// A preprocessor class derived from generic
// high-fidelity input generation class inputGen
// This derived class will implement methods to generate
// input files specific to EPIC 2016 code.
class NEMOSYS_EXPORT ep16Prep : public inputGen {
  public:
    explicit ep16Prep(const jsoncons::json &inputjson)
        : inputGen(inputjson), _shortForm(false), _mdb(nullptr) {};

    ~ep16Prep() override;

    static ep16Prep *readJSON(const std::string &jsonFName);
    static ep16Prep *readJSON(const jsoncons::json &inputjson);

  // basic methods
  public:
    void readJSON();
    void process() override;
    void read(const std::string &fname) override;
    void close(const std::string &fname) override;

  // feature specifications
  public:
    void addNdeSet() override {};
    void addElmSet() override {};

  // attribute assignment
  public:
    void addMat() override {};
    void addBC() override {};
    void addIC() override {};
    bool addMisc() override { return false; };

  // internal processing
  public:
    void wrtPre(const std::string &_tsk, const std::string &__tsk) override;
    void wrtCmnt(const std::string &cmnt) override;
    void wrtMsh(const std::string &_tsk, const std::string &__tsk) override;
    void wrtNdeCrds() override {};
    void wrtElmCon() override {};
    void wrtPost() override {};
    void wrtNdeSet() override {};
    void wrtElmSet() override {};
    void wrtMisc(const std::string &_tsk, const std::string &__tsk) override;
    void edit(const std::string &_tsk, const std::string &__tsk) override;

  private:
    std::string _fname;
    bool _shortForm;
    NEM::MSH::EXOMesh::exoMesh *_mdb;
    std::map<std::string, int> _mat;
    std::map<std::string, std::string> _bcs;
    std::stringstream _buffer;
};

#endif

#endif  // NEMOSYS_EP16PREP_H_
