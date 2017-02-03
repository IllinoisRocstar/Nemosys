# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "rbf_interp_nd.hpp"
# include "r8lib.hpp"
# include "rbfInterp.H"

int main ( );
void rbf_class_test01 ( );


//****************************************************************************80

int main ( )
{
  timestamp ( );
  cout << "\n";
  cout << "RBF Interpolation Class test:\n";
  cout << "  Test the RBFInterpolation library.\n";

  rbf_class_test01 ( );
//
//  Terminate.
//
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}

//****************************************************************************80

void rbf_class_test01 ( )

//****************************************************************************80
{
  double a;
  double app_error;
  double b;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  cout << "\n";
  cout << "RBF_Class_TEST01:\n";
  cout << "  RBF_WEIGHT computes weights for RBF interpolation.\n";
  cout << "  RBF_INTERP_ND evaluates the RBF interpolant.\n";
  cout << "  Use the multiquadratic basis function PHI1(R).\n";

  a = 0.0;
  b = 2.0;

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = new double[m*nd];

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  r8mat_transpose_print ( m, nd, xd, "  The product points:" );

  r0 = ( b - a ) / ( double ) ( n1d );

  cout << "\n";
  cout << "  Scale factor R0 = " << r0 << "\n";

  fd = new double[nd];
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  r8vec_print ( nd, fd, "  Function data:" );

  //w = rbf_weight ( m, nd, xd, r0, phi1, fd );
  //r8vec_print ( nd, w, "  Weight vector:" );

  // instantiate RBF Interpolant class
  RBFInterpolant* int1 = new RBFInterpolant(m, nd, GAUSS, 0.4);
  int1->setPointCoords(xd);
  int1->setPointData(fd);
//
//  #1: Interpolation test.  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  //fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );
  cout << "Starting the interpolation.\n";
  fi = int1->interpolate(ni, xi);
  r8vec_print ( ni, fi, "  Interpolated Function values:" );
  cout << "Finishing the interpolation.\n";
 

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] fi;
  delete [] xi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, seed );

  //fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );
  fi = int1->interpolate(ni, xi);

  fe = new double[ni];
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per 1000 samples = " << app_error << "\n";

  delete [] fd;
  delete [] fe;
  delete [] fi;
  delete [] x1d;
  delete [] xd;
  delete [] xi;

  return;
}
//****************************************************************************80
