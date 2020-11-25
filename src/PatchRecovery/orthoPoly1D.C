#include "orthoPoly1D.H"

#include <cmath>

// constructor
orthoPoly1D::orthoPoly1D(int _order, const std::vector<double> &x)
    : order(_order)
{
  a.resize(_order + 1);
  b.resize(_order);
  phi.resize(x.size(), _order + 1);
  phiTphiInv.resize(_order + 1, _order + 1);
  ComputePhiTPhiInv(x);
}

orthoPoly1D::orthoPoly1D(const orthoPoly1D &op)
{
  order = op.order;
  a = op.a;
  b = op.b;
  phi = op.phi;
  phiTphiInv = op.phiTphiInv;
}

orthoPoly1D &orthoPoly1D::operator=(const orthoPoly1D &op)
{
  if (this != &op)
  {
    order = op.order;
    a = op.a;
    b = op.b;
    phi = op.phi;
    phiTphiInv = op.phiTphiInv;
  }
  return *this;
}

double orthoPoly1D::EvaluateOrthogonal(int power, double xk) const
{
  double p0 = 1.;
  if (power == 0)
    return p0;
  double p1 = xk - a[0];
  if (power == 1)
    return p1;
  double p2 = std::numeric_limits<double>::quiet_NaN();
  for (int i = 1; i < power; ++i)
  {
    p2 = (xk - a[i]) * p1 - b[i - 1] * p0;
    p0 = p1;
    p1 = p2;
  }
  return p2;
}

/* Compute coefficients a_m, b_m to be used in recurrence relation for
   calculating orthogonal polynomial:
   a_m+1 = (sum_k=0^n-1 x_k*p_m^2(x_k))/(sum_k=0^n-1 p_m^2(x_k)) 
   b_m = (sum_k=0^n-1 p_m^2(x_k))/(sum_k=0^n-1 p_m-1^2(x_k)) */
void orthoPoly1D::ComputeAB(const std::vector<double> &x)
{
  a[0] = 0.;
  int n = x.size();
  for (int i = 0; i < n; ++i)
  {
    a[0] += x[i];
  }
  a[0] /= n;

  for (int i = 1; i < order; ++i)
  {
    double sum0, sum1, sum2;
    sum0 = sum1 = sum2 = 0.0;
    for (int j = 0; j < n; ++j)
    {
      double tmp0 = EvaluateOrthogonal(i - 1, x[j]);
      double tmp1 = EvaluateOrthogonal(i, x[j]);
      sum0 += tmp0 * tmp0;
      sum1 += tmp1 * tmp1;
      sum2 += x[j] * tmp1 * tmp1;
    }
    a[i] = sum2 / sum1;
    b[i - 1] = sum1 / sum0;
  }
}

/* Evaluate Orthogonal polynomials at data x and collect them in basis
   matrix phi */
void orthoPoly1D::EvaluateOrthogonals(const std::vector<double> &x)
{
  ComputeAB(x);
  for (int i = 0; i < order + 1; ++i)
  {
    for (int j = 0; j < x.size(); ++j)
      phi(j, i) = EvaluateOrthogonal(i, x[j]);
  }
}

/* Compute inverted matrix for use in normal equation */
void orthoPoly1D::ComputePhiTPhiInv(const std::vector<double> &x)
{
  EvaluateOrthogonals(x);
  VectorXd tmp = (phi.transpose() * phi).diagonal();
  for (int i = 0; i < tmp.rows(); ++i)
    tmp(i) = 1. / tmp(i);
  phiTphiInv = tmp.asDiagonal();
}
