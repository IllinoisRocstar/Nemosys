#ifndef vtkMesh_H
#define vtkMesh_H

#include <meshBase.H>

using namespace nemAux;

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
    std::vector<std::vector<double>> getCellVec(int id);
    // get diameter of circumsphere of each cell
    std::vector<double> getCellLengths();
    // get center of a cell
    std::vector<double> getCellCenter(int cellID);
 
  // diagnostics
  public:
		void report();
		void write();
		void write(std::string fname); 

  // set and get point and cell data
  public:
    // set point data (numComponets per point determined by dim of data[0] 
    void setPointDataArray(const char* name, const std::vector<std::vector<double>>& data);
    // set cell data (numComponents per cell determined by dim of data[0])
    void setCellDataArray(const char* name, const std::vector<std::vector<double>>& data);
    // set scalar cell data
    void setCellDataArray(const char* name, const std::vector<double>& data);
    // remove point data with given id from dataSet if it exists
    void unsetPointDataArray(int arrayID);
    void unsetPointDataArray(const char* name);
    // remove cell data with given id from dataSet if it exists
    void unsetCellDataArray(int arrayID);
    void unsetCellDataArray(const char* name);
    // remove field data with given id from dataSet
    void unsetFieldDataArray(const char* name);

};



#endif
