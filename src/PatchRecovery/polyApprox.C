#include "PatchRecovery/polyApprox.H"

#include <iostream>

#include <Eigen/LU>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector4d;
typedef Eigen::Matrix<double, 10, 1> Vector10d;


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
