#ifndef meshPhys_H
#define meshPhys_H

// Nemosys
#include <vtkAnalyzer.H>

//MAd Lib
#include <MAdLib.h>
#include <ModelInterface.h>
#include <MeshDataBaseInterface.h>
#include <GmshEntities.h>

// stl
#include <functional>

/* holds relevant info for point data and provides access interface*/
class PointDataArray
{
  public:
  // constructors
    PointDataArray( std::string _name, int _numComponent, int _numTuple, 
                    const std::vector<std::vector<double>>& _pntData);

    ~PointDataArray(){};

  public:
  // operators, access
  	double& operator()(int i, int j) { return pntData[i*numComponent +j ]; }
  	const double& operator()(int i, int j) const { return pntData[i*numComponent + j]; }  
  	int getNumComponent() { return numComponent; }
  	int getNumTuple() { return numTuple; }
  	std::string getName() { return name; }
    // get 1d pntData ({p11,p12,...,pn1,pn2,...})
    std::vector<double> getData() { return pntData; }
    // get folded pntData ({{p11,p12,...},...,{pn1,pn2,...})
    std::vector<std::vector<double>> getFoldData();  

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
    // constructs instance with vtkAnalyzer instance and populates pntData
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
  //access stuff

    // get array i from pntData
    PointDataArray& getPointData(int i) { return pntData[i]; }
    // get numComponents of array i from pntData
    int getDimArray(int i) { return pntData[i].getNumComponent(); }

  public:
  // relevant computations (gradients etc.)

    // compute gradient of point data at cell
    std::vector<double> ComputeGradAtCell(int cell, int array);
    // compute value of point data at center of cell
    std::vector<double> ComputeValAtCell(int cell, int array);
    // compute L2 norm of gradient of point data at each cell 
    std::vector<double> ComputeL2GradAtAllCells(int array);
    // compute value of point data at center of each cell
    std::vector<double> ComputeValAtAllCells(int array);
    // compute L2 norm of value of point data at center of each cell
    std::vector<double> ComputeL2ValAtAllCells(int array);
    // get diameter of circumsphere of each cell
    std::vector<double> GetCellLengths();

  
  public:
  // write functions and mesh processing/adaptation
    
    // writes a background mesh with sizes defined by 
    // point data interpolated to cell center
    void writeBackgroundMSH(string filename, std::vector<double> sizes);

    // generate background size field based on values or gradient
    // and write it to msh file named backgroundSF.msh
    void createSizeField(int array_id, std::string method, 
                                        double dev_mult);

    /* This function essentially does the following:

    1) Registers the relevant data with the nodal data manager
    2) Runs the adapter
    3) Unclassifies the boundary elements for proper output
    4) Writes the refined mesh to a .msh file without data
    5) Gets the data after refinement
    6) Converts the refined mesh to a vtk file without data
    7) Writes the post-refinement data to a refined vtk and gmsh mesh 
    8) Misc: Memory management, string trimming/file naming      */
    void Refine(MAd::MeshAdapter* adapter, MAd::pMesh& mesh,
                int array_id, int dim, std::string outMeshFile);    

  private:
    std::vector<PointDataArray> pntData ; // all pointData on mesh
      
};

//---------------------Auxiliarry Functions---------------------------------//

// flattens vector of vectors
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v);

// folds vector into vector of vectors
template <typename T>
std::vector<std::vector<T>> fold(const std::vector<T>& v, int dim);

// adds two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& x, 
                         const std::vector<T>& y);

// compute L2 norm of vec
double L2_Norm(const std::vector<double>& x);

// computes reciprocal of argument
template <typename T>
T reciprocal(T x);

// computes reciprocal of vector elements, in place
template <typename T>
void reciprocal_vec(std::vector<T>& x);

// find minmax of vector, excluding inf 
std::vector<double> getMinMax(const std::vector<double>& x);

// scales x from range [xmin, xmax] to within range [ymin, ymax]
// if x is inf, the scaled value will be ymax
double scale_to_range(double x, const std::vector<double>& xminmax, 
                      const std::vector<double>& yminmax);

// scales each number in vector in range [xmin, xmax] to [ymin,ymax]
void scale_vec_to_range(std::vector<double>& x, 
                        const std::vector<double>& xminmax,
                        const std::vector<double>& yminmax);

// multiplies vector by scalar
template <typename T>
std::vector<T> operator*(T a, std::vector<T>& x);

// get average and stdev of values
std::vector<double> getMeanStdev(std::vector<double>& x);

// generates boolean array with 1 if value >= tol, 0 otherwise
std::vector<int> cellsToRefine(std::vector<double>& values, double tol);

// string trimming for consistent file names 
std::string trim_fname(std::string name, std::string ext);
//-------------------------------------------------------------------------//

#endif
