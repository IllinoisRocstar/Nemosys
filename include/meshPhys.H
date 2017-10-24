#ifndef meshPhys_H
#define meshPhys_H

#include <vtkAnalyzer.H>
#include <MAdLib.h>

// stl

// Auxiliarry Functions
#include <functional>

// flattens vector of vectors
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v)
{
	std::size_t size = 0;
	for (const auto& sub : v)
		size+= sub.size();
	std::vector<T> result;
	result.reserve(size);
	for (const auto& sub : v)
		result.insert(result.end(), sub.begin(), sub.end());
	return result;
}

// adds two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& x, 
                              const std::vector<T>& y)
{
  if (x.size() != y.size())
  {
    std::cout << "vectors must be same length for addition" << std::endl;
    exit(1);
  }
  
  std::vector<T> result;
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), y.begin(), 
                 std::back_inserter(result), std::plus<T>());
  return result;
} 

// multiplies vector by scalar, in place
template <typename T>
void operator*(T a, std::vector<T>& x)
{
  std::transform(x.begin(), x.end(), x.begin(), std::bind1st(std::multiplies<T>(),a));
}  

// computes reciprocal of argument
template <typename T>
T reciprocal(T x)
{
  return (T) 1/x ;
}

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x)
{
  std::transform(x.begin(), x.end(), x.begin(), reciprocal<T>);
}

// compute L2 norm of vec
double L2_Norm(const std::vector<double>& x)
{
  double result = 0.0;
  for (int i = 0; i < x.size(); ++i)
    result += x[i]*x[i];
  
  return std::sqrt(result);
}

std::vector<double> getMinMax(const std::vector<double>& x)
{
  std::vector<double> result(2);
  auto minmax = std::minmax_element(x.begin(), x.end());
  result[0] = *minmax.first;
  result[1] = *minmax.second;
  return result;
} 

// scales x from range [xmin, xmax] to within range [ymin, ymax]
double scale_to_range(double x, const std::vector<double>& xminmax, 
                                const std::vector<double>& yminmax)
{
  return yminmax[0] + (yminmax[1] - yminmax[0])*(x - xminmax[0])/(xminmax[1] - xminmax[0]);
}

void scale_vec_to_range(std::vector<double>& x, 
                        const std::vector<double>& xminmax,
                        const std::vector<double>& yminmax)
{
  for (int i = 0; i < x.size(); ++i)
    x[i] = scale_to_range(x[i], xminmax, yminmax);
}


/* holds relevant info for point data and provides access interface*/
class PointDataArray
{
  public:
  // constructors
    PointDataArray( std::string _name, int _numComponent, int _numTuple, 
                    const std::vector<std::vector<double>>& _pntData):
                    name(_name), numComponent(_numComponent), numTuple(_numTuple)
                    { pntData = flatten (_pntData) ; };

    ~PointDataArray(){};

  public:
  // operators, access
  	double& operator()(int i, int j) { return pntData[i*numComponent +j ]; }
  	const double& operator()(int i, int j) const { return pntData[i*numComponent + j]; }  
  	int getNumComponent() { return numComponent; }
  	int getNumTuple() { return numTuple; }
  	std::string getName() { return name; }
  

  private: 
    int numComponent, numTuple; // dim of data point, number of data points 
    std::string name;           // name of data array
    std::vector<double> pntData; // data array
};

/* meshPhys inherits from vtkAnalyzer and contains methods for
   evaluating physical data on the mesh. Its core private member
   is a vector of all pointData arrays on the mesh */
class meshPhys: public vtkAnalyzer
{
  public:
  // constructors
    meshPhys(char* fname):vtkAnalyzer(fname)
    {
      read();
      
      if (dataSet)
      {
        getNumberOfPointData(); // number of point data arrays
        getNumberOfPoints();    // number of points
        getNumberOfCells();     // number of cells
      
        if (numberOfPointData < 1)
        {
          std::cout << "No point data found!" << std::endl;
          exit(1);
        }

        for (int i = 0; i < numberOfPointData; ++i)  
        {
          int numTuple, numComponent;
          std::vector<std::vector<double>> pntData_tmp;
          getPointDataArray(i, pntData_tmp, numTuple, numComponent);
          std::string name = dataSet->GetPointData()->GetArrayName(i);
          pntData.push_back(PointDataArray(name, numComponent, numTuple, pntData_tmp));
        }
      } 
    };

    ~meshPhys(){};

  public:
  //access to pntData
    PointDataArray& getPointData(int i) { return pntData[i]; }

  public:
  // relevant computations (gradients etc.)

  std::vector<double> ComputeGradAtCell(int cell, int array);
  std::vector<double> ComputeL2GradAtAllCells(int array);
  std::vector<double> ComputeGradAtPoint(int pnt, int array);
  std::vector<double> GetCellLengths();

  
  public:
  // write functions
    
    void writeBackgroundMSH(string filename, std::vector<double> sizes);


  private:
    std::vector<PointDataArray> pntData ;
      
};

#endif
