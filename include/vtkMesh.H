#ifndef vtkMesh_H
#define vtkMesh_H

#include <meshBase.H>

class vtkMesh : public meshBase
{
  // constructor and destructor
  public:

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
     
  // diagnostics
  public:
    void report(char* fname);
    void write(std::string fname); 
  
  // set point and cell data
  public:

    void setPointDataArray(const char* name, std::vector<std::vector<double>>& data);
    void setCellDataArray(const char* name, std::vector<std::vector<double>>& data);



};



#endif
