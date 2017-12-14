#ifndef MESHQUALITY_H
#define MESHQUALITY_H

#include <meshBase.H>

class MeshQuality
{
  public:
    MeshQuality(){}
    MeshQuality(meshBase* _mesh);
    ~MeshQuality();

  public:
    void checkMesh(std::ostream& outputStream);
    void checkMesh();
    void checkMesh(std::string fname); vtkSmartPointer<vtkDoubleArray> getStats(int n);

  private:
    meshBase* mesh;
    vtkSmartPointer<vtkMeshQuality> qualityFilter;

};


#endif
