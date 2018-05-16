#ifndef MESHGEN_H
#define MESHGEN_H

#include <string>
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>

// method to write vtk grids
template<class TWriter> void writeVTFile(const std::string& fname, 
                                         vtkSmartPointer<vtkDataSet> dataSet)
{
  vtkSmartPointer<TWriter> Writer = vtkSmartPointer<TWriter>::New();
  Writer->SetFileName(&fname[0u]);
  Writer->SetInputData(dataSet);
  Writer->Write();
}


class meshingParams
{
  public:
    meshingParams(){};
    virtual ~meshingParams(){};    
    
};

class meshGen 
{
  public:
    meshGen():dataSet(0){}
    virtual ~meshGen(){}
    
    // creates generator with default parameters
    static meshGen* Create(std::string fname, std::string meshEngine);
    // creates generater with specified parameters
    static meshGen* Create(std::string fname, std::string meshEngine, meshingParams* params);
    virtual int createMeshFromSTL(const char* fname) = 0;
    vtkSmartPointer<vtkDataSet> getDataSet(); 
 
  protected:
    vtkSmartPointer<vtkDataSet> dataSet;    

};


#endif
