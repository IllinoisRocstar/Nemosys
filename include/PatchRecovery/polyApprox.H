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
#ifndef NEMOSYS_POLYAPPROX_H_
#define NEMOSYS_POLYAPPROX_H_

#include "nemosys_export.h"

#include <vector>
#include <memory>

#include <Eigen/Core>

class NEMOSYS_EXPORT polyApprox
{
  public:
    // assumes coordinates are preconditioned
    polyApprox(int _order,
               const std::vector<std::vector<double>> &coords);

    ~polyApprox() = default;

    // disable copy and copy-assignment constructors
    polyApprox(const polyApprox &) = delete;
    polyApprox &operator=(const polyApprox &) = delete;

    // compute coefficients for polynomial expansion of sampled function
    void computeCoeff(const Eigen::VectorXd &data);
    void resetCoeff();
    double eval(const std::vector<double> &coord) const;

    static std::unique_ptr<polyApprox>
    CreateUnique(int order,
                 const std::vector<std::vector<double>> &coords);

  private:
    int order;
    // matrix of basis polynomials evaluated at all coords
    std::vector<Eigen::VectorXd> basis;
    // lhs matrix in normal equation (phiTphi*a = phiT*f <=> A*a = b)
    Eigen::MatrixXd A;
    // rhs vector in normal equation
    Eigen::VectorXd b;
    // coefficients of polynomial approximant
    Eigen::VectorXd a;

    Eigen::VectorXd computeBasis(const std::vector<double> &coord) const;
};

#endif  // NEMOSYS_POLYAPPROX_H_
