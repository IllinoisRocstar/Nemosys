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

};



#endif
