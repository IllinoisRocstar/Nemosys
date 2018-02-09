#ifndef POLYAPPROX_H
#define POLYAPPROX_H

#include <iostream>
#include <ostream>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
class polyApprox
{
  public:
    // assumes coordinates are preconditioned
    polyApprox(const int _order, const std::vector<std::vector<double>>& _coords);
    ~polyApprox(){}
    
    // compute coefficients for polynomial expansion of sampled function
    void computeCoeff(const VectorXd data);
    double eval(const std::vector<double>& coord);
  
  private:
    std::vector<std::vector<double>> coords;
    VectorXd data;
    int order;
    // matrix of basis polynomials evaluated at all coords
    std::vector<VectorXd> basis;
    // lhs matrix in normal equation (phiTphi*a = phiT*f <=> A*a = b)
    MatrixXd A;
    // lhs vector in normal equation
    VectorXd b;
    // coefficients of polynomial approximant
    VectorXd a;
    
    VectorXd computeBasis(const std::vector<double>& coord);


    // private lvalue copy ctor disables it (not implemented)
		polyApprox(const polyApprox&) = delete;
    // private lvalue copy assignment ctor disables it (not implemented)
    polyApprox& operator=(const polyApprox&) = delete; 
    
    
};


#endif
