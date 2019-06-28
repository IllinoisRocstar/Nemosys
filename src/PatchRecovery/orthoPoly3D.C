#include "orthoPoly3D.H"

#ifdef NEMOSYS_DEBUG
  #include <AuxiliaryFunctions.H>
#endif

// remove given column from matrix in place (VERY SLOW)
void removeColumn(MatrixXd &matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
        matrix.rightCols(numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}


// remove column from matrix, out of place
// is efficient given referenced matrices are stored as column-major
void removeColumn(const MatrixXd &matrix,
                  MatrixXd &matrix_red,
                  const std::vector<unsigned int> &toRemove)
{
  int m = 0;
  for (int i = 0; i < matrix.cols(); ++i)
  {
    bool remove = false;
    for (unsigned int k : toRemove)
    {
      if (i == k)
      {
        remove = true;
        break;
      }
    }
    if (!remove)
    {
      matrix_red.col(m) = matrix.col(i);
      m += 1;
    }
  }
}


// exclude columns from a matrix and construct the transpose
// of such an exclusion. This method is fast given the immutable
// matrix is stored as column-major and the mutable one is row-major
void removeRowT(const MatrixXd &matrix,
                MatrixXdRM &matrix_red,
                const std::vector<unsigned int> &toRemove)
{
  int m = 0;
  for (int i = 0; i < matrix.cols(); ++i)
  {
    bool remove = false;
    for (unsigned int k : toRemove)
    {
      if (i == k)
      {
        remove = true;
        break;
      }
    }
    if (!remove)
    {
      matrix_red.row(m) = matrix.col(i);
      m += 1;
    }
  }
}

// default ctor
orthoPoly3D::orthoPoly3D() : finished(false),
                             opx(new orthoPoly1D()),
                             opy(new orthoPoly1D()),
                             opz(new orthoPoly1D())
{
  a.resize(0);
}

void orthoPoly3D::initCheck()
{
  if (order == 1)
    toRemove = {3, 5, 6, 7};
  else if (order == 2)
    toRemove = {5, 7, 8, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25,
                26};
  else if (order == 3)
    toRemove = {7, 10, 11, 13, 14, 15, 19, 22, 23, 25, 26,
                27, 28, 29, 30, 31, 34, 35, 37, 38, 39, 40,
                41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52,
                53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
  else
  {
    std::cerr << "Polynomial order greater than 3 is not supported"
              << std::endl;
    exit(1);
  }
}

// complete ctor
orthoPoly3D::orthoPoly3D(int _order, const VectorXd &sigma,
                         const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &z)
    : order(_order), finished(false)
{
  initCheck();

  opx = std::unique_ptr<orthoPoly1D>(new orthoPoly1D(order, x));
  opy = std::unique_ptr<orthoPoly1D>(new orthoPoly1D(order, y));
  opz = std::unique_ptr<orthoPoly1D>(new orthoPoly1D(order, z));
  computeA(sigma);
}

orthoPoly3D::orthoPoly3D(int _order,
                         const std::vector<std::vector<double>> &coords)
    : finished(false), order(_order)
{
  //int root(round(cbrt(coords.size())));
  //if (coords.size() != root * root * root &&
  //    coords.size() != (root + 1) * (root + 1) * (root + 1))
  //{
  //  std::cerr << "coords are definitely not in a cartesian product set"
  //            << std::endl;
  //  exit(1);
  //}

  initCheck();

  std::vector<double> x(coords.size());
  std::vector<double> y(coords.size());
  std::vector<double> z(coords.size());

  for (int i = 0; i < coords.size(); ++i)
  {
    x[i] = coords[i][0];
    y[i] = coords[i][1];
    z[i] = coords[i][2];
  }
  std::cout << "non-unique coords" << std::endl;
#ifdef NEMOSYS_DEBUG
  nemAux::printVec(x); //nemAux::printVec(y); nemAux::printVec(z);
#endif
  //std::vector<double> x_unique;
  //for (int i = 0; i < x.size(); ++i)
  //{
  //  x_unique.push_back(x[i]);
  //  for (int j = i + 1; j < x.size(); ++j)
  //  {
  //    if (x_unique[i] != x[j])
  //  }
  //}

  //std::sort(x.begin(), x.end());
  //x.erase(std::unique(x.begin(), x.end(), Pred), x.end());
  //std::cout << "unique coords" << std::endl;
  //nemAux::printVec(x); //nemAux::printVec(y); nemAux::printVec(z);

  opx = std::unique_ptr<orthoPoly1D>(new orthoPoly1D(order, x));
  opy = std::unique_ptr<orthoPoly1D>(new orthoPoly1D(order, y));
  opz = std::unique_ptr<orthoPoly1D>(new orthoPoly1D(order, z));
}

// move ctor
orthoPoly3D::orthoPoly3D(orthoPoly3D &&op)
    : finished(op.finished),
      opx(std::move(op.opx)),
      opy(std::move(op.opy)),
      opz(std::move(op.opz)),
      a(op.a),
      toRemove(op.toRemove),
      order(op.order)
{
  op.Reset();
}


// move assignment  
orthoPoly3D &orthoPoly3D::operator=(orthoPoly3D &&op)
{
  // check against self assignment
  if (this != &op)
  {
    // release held resources
    this->Reset();
    // move resources from op to this (i.e. pilfer)
    opx = std::move(op.opx);
    opy = std::move(op.opy);
    opz = std::move(op.opz);
    a = op.a;
    toRemove = op.toRemove;
    finished = op.finished;
    order = op.order;
    // set op to default state  
    op.Reset();
  }
  return *this;
}

orthoPoly3D *
orthoPoly3D::Create(int _order, const std::vector<std::vector<double>> &coords)
{
  return new orthoPoly3D(_order, coords);
}

std::unique_ptr<orthoPoly3D>
orthoPoly3D::CreateUnique(int _order,
                          const std::vector<std::vector<double>> &coords)
{
  return std::unique_ptr<orthoPoly3D>(
      orthoPoly3D::Create(_order, coords));
}


// compute coefficients for polynomial expansion of sampled function
void orthoPoly3D::computeA(const VectorXd &sigma)
{
  if (!finished)
  {
#ifdef NEMOSYS_DEBUG
    nemAux::Timer T;
#endif

    MatrixXd tmp1 = opx->phi * opx->phiTphiInv;
    MatrixXd tmp2 = opy->phi * opy->phiTphiInv;
    MatrixXd tmp3 = opz->phi * opz->phiTphiInv;
#ifdef NEMOSYS_DEBUG
    nemAux::Timer T1;
    T1.start();
#endif
    MatrixXd kronProd_full = kroneckerProduct(tmp1, tmp2);
    kronProd_full = kroneckerProduct(kronProd_full, tmp3).eval();
#ifdef NEMOSYS_DEBUG
    T1.stop();
    std::cout << "Time for building kronprod in constructor: " << T1.elapsed()
              << "\n";
    std::cout << "kronProd_full: " << kronProd_full.rows() << " "
              << kronProd_full.cols() << "\n";
    std::cout << "Removing rows: ";
    nemAux::printVec(toRemove);
#endif
    MatrixXdRM kronProd_red(kronProd_full.cols() - toRemove.size(),
                            kronProd_full.rows());
#ifdef NEMOSYS_DEBUG
    T.start();
#endif
    removeRowT(kronProd_full, kronProd_red, toRemove);
    std::cout << "kronProd_red: " << kronProd_red.rows() << " "
              << kronProd_red.cols() << "\n";
    std::cout << "sigma size: " << sigma.size() << "\n";
    a = kronProd_red * sigma;
#ifdef NEMOSYS_DEBUG
    T.stop();
    std::cout << "Time: " << T.elapsed() << std::endl;
    std::cout << "shape of a: " << a.rows() << " " << a.cols() << std::endl;
#endif
    finished = true;
  }
}

void orthoPoly3D::resetA()
{
  if (finished)
  {
    a.resize(0);
    finished = false;
  }
}


const double orthoPoly3D::eval(const std::vector<double> &coord)
{
  if (!finished)
  {
    std::cerr << "Coefficients to basis have not been generated" << std::endl;
    exit(1);
  }

  MatrixXd phix(1, order + 1);
  MatrixXd phiy(1, order + 1);
  MatrixXd phiz(1, order + 1);
  for (int i = 0; i < order + 1; ++i)
  {
    phix(0, i) = opx->EvaluateOrthogonal(i, coord[0]);
    phiy(0, i) = opy->EvaluateOrthogonal(i, coord[1]);
    phiz(0, i) = opz->EvaluateOrthogonal(i, coord[2]);
  }

  MatrixXd Phi = kroneckerProduct(phix, phiy);
  Phi = kroneckerProduct(Phi, phiz).eval();
#ifdef NEMOSYS_DEBUG
  std::cout << Phi.rows() << " " << Phi.cols() << std::endl;
#endif
  MatrixXd Phi_red(Phi.rows(), Phi.cols() - toRemove.size());
  removeColumn(Phi, Phi_red, toRemove);

#ifdef NEMOSYS_DEBUG
  std::cout << Phi_red.rows() << " " << Phi_red.cols() << std::endl;
  std::cout << a.rows() << " " << a.cols() << std::endl;
#endif

  return (Phi_red * a)(0, 0);
}

// evaluate the fitting polynomial at coord
const double orthoPoly3D::operator()(const std::vector<double> &coord)
{
  return eval(coord);
}

#ifdef NEMOSYS_DEBUG
// compute multivariate polynomial approximation of prescribed order to sigma
void orthoPoly3D::run1(const VectorXd &sigma)
{
  MatrixXd NewPhi;
  NewPhi = kroneckerProduct(opx->phi, opy->phi);
  NewPhi = kroneckerProduct(NewPhi, opz->phi).eval();

  nemAux::Timer T;
  std::cout << "Removing stuff" << std::endl;
  T.start();

  MatrixXd NewPhi_red(NewPhi.rows(), NewPhi.cols() - toRemove.size());

  removeColumn(NewPhi, NewPhi_red, toRemove);

  T.stop();
  std::cout << "Time spent removing: " << T.elapsed() << std::endl;

  VectorXd sigma_approx = NewPhi_red * a;

  VectorXd error = (sigma_approx - sigma).cwiseAbs();
  std::cout << "max error " << error.maxCoeff() << std::endl;
  std::cout << "min error " << error.minCoeff() << std::endl;
}
#endif

void orthoPoly3D::Reset()
{
  if (opx) opx.reset();
  if (opy) opy.reset();
  if (opz) opz.reset();
  a.resize(0);
  finished = false;
  toRemove.resize(0);
}
