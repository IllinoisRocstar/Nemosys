#ifndef vtkMesh_H
#define vtkMesh_H

#include <meshBase.H>

class vtkMesh : public meshBase
{
  // constructor and destructor
  public:
    
    vtkMesh(){};
    vtkMesh(const char* fname);
    ~vtkMesh()
    {
      std::cout << "vtkMesh destroyed" << std::endl;
    } 

  // access
  public:
    // get point with id
    std::vector<double> getPoint(int id);
    // get cell with id : returns point indices and respective coordinates
    std::map<int, std::vector<double>> getCell(int id);
	  // returns cells sharing point (topological search)
	  std::vector<int> getCellsWithPoint(int pnt);
 
 
  // diagnostics
  public:
    void report(char* fname);
    void write(std::string fname, std::string extension); 
  
  // set and get point and cell data
  public:
    // set point data (numComponets per point determined by dim of data[0] 
    void setPointDataArray(const char* name, std::vector<std::vector<double>>& data);
    // set cell data (numComponents per cell determined by dim of data[0])
    void setCellDataArray(const char* name, std::vector<std::vector<double>>& data);
    // remove point data with given id from target if it exists
    void unsetPointDataArray(int arrayID);
    void unsetPointDataArray(const char* name);
    // remove cell data with given id from target if it exists
    void unsetCellDataArray(int arrayID);
    void unsetCellDataArray(const char* name);
    // transfer point data with given id from source to target


};



#endif
