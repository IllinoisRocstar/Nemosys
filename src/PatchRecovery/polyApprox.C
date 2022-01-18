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
#include "PatchRecovery/polyApprox.H"

#include <iostream>

#include <Eigen/LU>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector4d;
using Vector10d = Eigen::Matrix<double, 10, 1>;


polyApprox::polyApprox(const int _order,
                       const std::vector<std::vector<double>> &coords)
    : order(_order)//, coords(_coords)
{
  switch (order)
  {
    case 1:
    {
      A = MatrixXd::Zero(4, 4);
      b = VectorXd::Zero(4);
      a.resize(4);
      break;
    }
    case 2:
    {
      A = MatrixXd::Zero(10, 10);//.resize(10, 10);
      b = VectorXd::Zero(10);//resize(10);
      a.resize(10);
      break;
    }
    default:
    {
      std::cerr << "Error: order: " << order << " is not supported"
                << std::endl;
      exit(1);
    }
  }
  basis.resize(coords.size());
  for (int i = 0; i < coords.size(); ++i)
  {
    basis[i] = computeBasis(coords[i]);
    A = A + basis[i] * (basis[i].transpose()); // * basis[i];
  }
}


std::unique_ptr<polyApprox>
polyApprox::CreateUnique(const int order,
                         const std::vector<std::vector<double>> &coords)
{
  return std::unique_ptr<polyApprox>(new polyApprox(order, coords));
}


void polyApprox::computeCoeff(const VectorXd &data)
{
  for (int i = 0; i < basis.size(); ++i)
  {
    b = b + basis[i] * data(i); // basis[i].transpose() * data(i);
  }
  a = A.partialPivLu().solve(b);
}


void polyApprox::resetCoeff()
{
  a.setZero();
  b.setZero();
}


double polyApprox::eval(const std::vector<double> &coord) const
{
  VectorXd P = computeBasis(coord);
  return (P.transpose() * a)(0, 0);
}


VectorXd polyApprox::computeBasis(const std::vector<double> &coord) const
{
  switch (order)
  {
    case 1:
    {
      Vector4d basisVec;
      basisVec(0) = 1;
      basisVec(1) = coord[0];
      basisVec(2) = coord[1];
      basisVec(3) = coord[2];
      return basisVec;
    }
    case 2:
    {
      Vector10d basisVec;
      basisVec(0) = 1;
      basisVec(1) = coord[0];
      basisVec(2) = coord[1];
      basisVec(3) = coord[2];
      basisVec(4) = coord[0] * coord[0];
      basisVec(5) = coord[0] * coord[1];
      basisVec(6) = coord[0] * coord[2];
      basisVec(7) = coord[1] * coord[1];
      basisVec(8) = coord[1] * coord[2];
      basisVec(9) = coord[2] * coord[2];
      return basisVec;
    }
    default:
    {
      std::cerr << "Error: order " << order << " is not supported."
                << std::endl;
      exit(1);
    }
  }
}
