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
#ifndef NEMOSYS_ORTHOPOLY1D_H_
#define NEMOSYS_ORTHOPOLY1D_H_

#include "nemosys_export.h"

#include <iostream>
#include <vector>

#include <Eigen/Core>

class NEMOSYS_EXPORT orthoPoly1D
{
  // ctors and dtor
  public:
    orthoPoly1D(int _order, const std::vector<double> &x);
    orthoPoly1D() = default;

    orthoPoly1D(const orthoPoly1D &op);
    orthoPoly1D &operator=(const orthoPoly1D &op);

    ~orthoPoly1D() = default;

  // methods invoked in constructor
  public:
    /* Use recurrence relation to calculate power^th orthogonal polynomial:
       p_m+1(x) = (x-a_m+1(x))*p_m(x) - b_m*p_m-1(x) for m>=1;
       p_0(x) = 1; p_1(x) = x - a[0] */
    double EvaluateOrthogonal(int power, double xk) const;

    /* Compute coefficients a_m, b_m to be used in recurrence relation for
       calculating orthogonal polynomial:
       a_m+1 = (sum_k=0^n-1 x_k*p_m^2(x_k))/(sum_k=0^n-1 p_m^2(x_k)) 
       b_m = (sum_k=0^n-1 p_m^2(x_k))/(sum_k=0^n-1 p_m-1^2(x_k)) */
    void ComputeAB(const std::vector<double> &x);

    /* Evaluate Orthogonal polynomials at data x and collect them in basis
       matrix phi */
    void EvaluateOrthogonals(const std::vector<double> &x);

    /* Compute inverted matrix for use in normal equation */
    void ComputePhiTPhiInv(const std::vector<double> &x);

  public:
    // polynomial order
    int order;
    // coordinate array
    //const std::vector<double> &x;
    // coefficient array to calculate value of orthogonal polynomial 
    std::vector<double> a;
    // coefficient array to calculate value of orthogonal polynomial 
    std::vector<double> b;
    // Matrix of polynomial basis phi, evaluated at x
    Eigen::MatrixXd phi;
    // Inverted matrix in normal equations a = [(phiT*phi)^-1]*phiT*f
    Eigen::MatrixXd phiTphiInv;
};

#endif  // NEMOSYS_ORTHOPOLY1D_H_
