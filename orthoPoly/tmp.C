#include <iostream>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <chrono>
#include <vector>

//using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;
//using Eigen::ArrayXd;
//using Eigen::KroneckerProduct;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdRM;

class Timer {
private:
  typedef std::chrono::time_point<std::chrono::system_clock> time_t;

public:
  Timer() : startTime(), stopTime() {}

  time_t start()   { return (startTime = std::chrono::system_clock::now()); }
  time_t stop()    { return (stopTime  = std::chrono::system_clock::now()); }
  double elapsed() { return std::chrono::duration_cast<std::chrono::milliseconds>
                                                      (stopTime-startTime).count(); }

private:
  time_t startTime, stopTime;
};


// remove given row from matrix
void removeRow(MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = 
          matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

// remove given column from matrix
void removeColumn(MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = 
          matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void removeRow(const MatrixXd& matrix, 
               MatrixXd& matrix_red, std::vector<unsigned int>& toRemove)
{
  int numRows = matrix.rows()-toRemove.size();
  int numCols = matrix.cols();
  int m,n;
  for (int j = 0; j < matrix.cols(); ++j)
  {
    m = 0;
    for (int i = 0; i < matrix.rows(); ++i)
    {
      bool remove = 0;
      for (int k = 0; k < toRemove.size(); ++k)
      {
        if (i == toRemove[k])
        {
          remove = 1;
          break;
        }
      }
      if (!remove)
      {
        matrix_red(m,j) = matrix(i,j);
        m+=1;
      }
    }
  }
}

void removeRowNew(const MatrixXd& matrix,
                  MatrixXdRM& matrix_red, 
                  std::vector<unsigned int>& toRemove)
{

  int m = 0;
  for (int i = 0; i < matrix.rows(); ++i)
  {
    bool remove = 0;
    for (int k = 0; k < toRemove.size(); ++k)
    {
      if (i == toRemove[k])
      { 
        remove = 1;
        break;
      }
    }
    if (!remove)
    {
      matrix_red.row(m) = matrix.row(i);
      m+=1;
    } 
  }
}

void removeRowNew(const MatrixXdRM& matrix,
                  MatrixXdRM& matrix_red, 
                  std::vector<unsigned int>& toRemove)
{

  int m = 0;
  for (int i = 0; i < matrix.rows(); ++i)
  {
    bool remove = 0;
    for (int k = 0; k < toRemove.size(); ++k)
    {
      if (i == toRemove[k])
      { 
        remove = 1;
        break;
      }
    }
    if (!remove)
    {
      matrix_red.row(m) = matrix.row(i);
      m+=1;
    } 
  }
}
void removeRowNew(const MatrixXd& matrix,
                  MatrixXd& matrix_red, 
                  std::vector<unsigned int>& toRemove)
{

  int m = 0;
  for (int i = 0; i < matrix.rows(); ++i)
  {
    bool remove = 0;
    for (int k = 0; k < toRemove.size(); ++k)
    {
      if (i == toRemove[k])
      { 
        remove = 1;
        break;
      }
    }
    if (!remove)
    {
      matrix_red.row(m) = matrix.row(i);
      m+=1;
    } 
  }
}


MatrixXd removeRow(const MatrixXd& matrix, 
               std::vector<unsigned int>& toRemove)
{
  int numRows = matrix.rows()-toRemove.size();
  int numCols = matrix.cols();
  MatrixXd result(numRows,numCols);
  int m,n;
  for (int j = 0; j < matrix.cols(); ++j)
  {
    m = 0;
    for (int i = 0; i < matrix.rows(); ++i)
    {
      bool remove = 0;
      for (int k = 0; k < toRemove.size(); ++k)
      {
        if (i == toRemove[k])
        {
          remove = 1;
          break;
        }
      }
      if (!remove)
      {
        result(m,j) = matrix(i,j);
        m+=1;
      }
    }
  }
  return result;
}

void removeColumn(const MatrixXd& matrix,
                  MatrixXd& matrix_red,
                  std::vector<unsigned int>& toRemove)
{
  int m = 0;
  for (int i = 0; i < matrix.cols(); ++i)
  {
    bool remove = 0;
    for (int k = 0; k < toRemove.size(); ++k)
    {
      if (i == toRemove[k])
      { 
        remove = 1;
        break;
      }
    }
    if (!remove)
    {
      matrix_red.col(m) = matrix.col(i);
      m+=1;
    } 
  }
  
}

int main()
{
  MatrixXd mat = MatrixXd::Random(64,1000000);
  MatrixXdRM matRM = MatrixXdRM::Random(64,1000000); 
  std::vector<unsigned int> toRemove = {7,10,11,13,14,15,19,22,23,25,26,
                                        27,28,29,30,31,34,35,37,38,39,40,
                                        41,42,43,44,45,46,47,49,50,51,52,
                                        53,54,55,56,57,58,59,60,61,62,63};
  
  Timer T;

  MatrixXdRM mat_red(mat.rows()-toRemove.size(),mat.cols());
  MatrixXdRM mat_red1(mat.rows()-toRemove.size(),mat.cols());
  
  MatrixXd mat_redold(mat.rows()-toRemove.size(),mat.cols());
  T.start();
  removeRowNew(mat, mat_redold, toRemove);
  T.stop();
  std::cout << "T old: " << T.elapsed() << std::endl;
  
  std::cout << mat.rows() << " " << mat.cols() << std::endl;
  //std::cout << mat.col(7) << std::endl;
  T.start();
  removeRowNew(mat, mat_red, toRemove);
  T.stop();
  std::cout << mat_red.rows() << " " << mat_red.cols() << std::endl;
  std::cout << "T new: " << T.elapsed() << std::endl;

  T.start();
  removeRowNew(matRM,mat_red1,toRemove);
  T.stop(); 
  std::cout << mat_red1.rows() << " " << mat_red1.cols() << std::endl;
  std::cout << "T new: " << T.elapsed() << std::endl;

  //MatrixXd mat = MatrixXd::Random(10000,10000);
  //VectorXd vec = VectorXd::Random(10000);

  //MatrixXd result(10000,1); 
  //Timer T;
  //T.start();
  //result = mat*vec;
  //T.stop();
  //std::cout << T.elapsed()*.001 << std::endl;
  return 0;


}
    

