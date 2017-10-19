#ifndef meshPhys_H
#define meshPhys_H

#include <vtkAnalyzer.H>
#include <MAdLib.h>

// stl


/* holds relevant info for point data and provides access interface*/
class PointDataArray
{
  public:
  // constructors
    PointDataArray( std::string _name, int _numComponent, int _numTuple, 
                    std::vector<std::vector<double>> _pntData):
                    name(_name), numComponent(_numComponent), numTuple(_numTuple)
                    { pntData = _pntData; };

    ~PointDataArray(){};

  public:
  // operators, access
  double& operator()(int i, int j) { return pntData[i][j]; }
  const double& operator()(int i, int j) const {return pntData[i][j]; }  
  int getNumComponent() { return numComponent; }
  int getNumTuple() { return numTuple; }
  std::string getName() { return name; }
  

  private: 
    int numComponent, numTuple; // dim of data point, number of data points 
    std::string name;           // name of data array
    std::vector<std::vector<double>> pntData; // data array
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

  /*  Algorithm: ComputeGradAtCenter:
      BEGIN:
        for each CELL in CELLS:
          JacboianInverse(double** inverse, double derivs[12]);

          |Now, we have the inverse jacobian (derivatives of shape function in iso-param space)
          |We use the interpolation derivatives (in iso-param space) to get the value of 
          |the functions derivative at cell center (in iso-param space)
          |Then, we use the inverse jacobian to transform back to the original coordinates

          |derivs[12] -> derivs[4][3] (deriv of shape function for each point)
          |Consider point data U=(ux,uy,uz)
          
          dU/dx = dU/dy = dU/dz = 0;
          for each POINT in CELL:
            dU/dx += derivs[POINT][0]*ux[POINT]; 
            dU/dy += derivs[POINT][1]*ux[POINT];
            dU/dy += derivs[POINT][2]*ux[POINT]; 

          [dU/dx dU/dy dU/dz]_orig = inverse*[dU/dx dU/dy dU/dz];              
          return [dU/dx dU/dy dU/dz]_orig; 
      END
  */

  std::vector<std::vector<double>> ComputeGradAtPoint(int cell, int array);


  // returns vector of parametric coords for all cells
  std::vector<double*> getParametricCoords();

  
 
  public:
  // write functions
    
    void writeBackgroundMSH(string filename, std::vector<double> sizes);


  private:
    std::vector<PointDataArray> pntData ;
      
};

#endif
