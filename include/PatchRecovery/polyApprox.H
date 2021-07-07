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
