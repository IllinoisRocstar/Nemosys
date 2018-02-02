#ifndef ORTHOPOLY1D_H
#define ORTHOPOLY1D_H

#include <iostream>
#include <ostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class orthoPoly1D
{
  // ctors and dtor
  public:
    orthoPoly1D(int _order, const std::vector<double>& x);
		orthoPoly1D(){}
		orthoPoly1D(const orthoPoly1D& op);
		orthoPoly1D& operator=(const orthoPoly1D& op);		
		~orthoPoly1D(){}

  // methods invoked in constructor
  public:
    
    /* Use recurrence relation to calculate power^th orthogonal polynomial:
       p_m+1(x) = (x-a_m+1(x))*p_m(x) - b_m*p_m-1(x) for m>=1;
       p_0(x) = 1; p_1(x) = x - a[0] */
    double EvaluateOrthogonal(int power, double xk);

    /* Compute coeficcients a_m,b_m to be used in recurrence relation for 
       calculating orthogonal polynomial:
       a_m+1 = (sum_k=0^n-1 x_k*p_m^2(x_k))/(sum_k=0^n-1 p_m^2(x_k)) 
       b_m = (sum_k=0^n-1 p_m^2(x_k))/(sum_k=0^n-1 p_m-1^2(x_k)) */
    void ComputeAB(const std::vector<double>& x);

    /* Evaluate Orthogonal polynomials at data x and collect them in basis
       matrix phi */
    void EvaluateOrthogonals(const std::vector<double>& x);
 
    /* Compute inverted matrix for use in normal equation */
    void ComputePhiTPhiInv(const std::vector<double>& x);


  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // polynomial order
    int order;
    // coordinate array
    //const std::vector<double>& x;
    // coefficient array to calculate value of orthogonal polynomial 
    std::vector<double> a;
    // coefficient array to calculate value of orthogonal polynomial 
    std::vector<double> b;
    // Matrix of polynomial basis phi, evaluated at x
    MatrixXd phi;
    // Inverted matrix in normal equations a = [(phiT*phi)^-1]*phiT*f
    MatrixXd phiTphiInv;
};
#endif
