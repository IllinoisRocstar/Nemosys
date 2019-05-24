#ifndef MESHQUALITY_H
#define MESHQUALITY_H

#include <meshBase.H>
#include <vtkMeshQuality.h>
#include <vtkDoubleArray.h>

class MeshQuality
{
  public:
    MeshQuality(){}
    MeshQuality(const meshBase* _mesh);
    ~MeshQuality();

  public:
    void checkMesh(std::ostream& outputStream);
    void checkMesh();
    void checkMesh(const std::string &fname);
    vtkSmartPointer<vtkDoubleArray> getStats(int n);

  private:
    const meshBase* mesh;
    vtkSmartPointer<vtkMeshQuality> qualityFilter;
  
    MeshQuality(const MeshQuality& that) = delete;
    MeshQuality& operator=(const MeshQuality& that) = delete; 
};


#endif
