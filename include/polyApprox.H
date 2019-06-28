#ifndef POLYAPPROX_H
#define POLYAPPROX_H

#include "nemosys_export.h"

#include <vector>
#include <memory>

// Eigen
#include <Eigen/Dense>


using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    void computeCoeff(const VectorXd &data);
    void resetCoeff();
    double eval(const std::vector<double> &coord) const;

    static std::unique_ptr<polyApprox>
    CreateUnique(int order,
                 const std::vector<std::vector<double>> &coords);

  private:
    int order;
    // matrix of basis polynomials evaluated at all coords
    std::vector<VectorXd> basis;
    // lhs matrix in normal equation (phiTphi*a = phiT*f <=> A*a = b)
    MatrixXd A;
    // rhs vector in normal equation
    VectorXd b;
    // coefficients of polynomial approximant
    VectorXd a;

    VectorXd computeBasis(const std::vector<double> &coord) const;
};

#endif
