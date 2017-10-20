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
  std::vector<double> ComputeGradAtPoint(int pnt, int array);

  
  public:
  // write functions
    
    void writeBackgroundMSH(string filename, std::vector<double> sizes);


  private:
    std::vector<PointDataArray> pntData ;
      
};

#endif
