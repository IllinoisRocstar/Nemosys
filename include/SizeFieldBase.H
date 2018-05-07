#ifndef SIZEFIELDGEN_H
#define SIZEFIELDGEN_H

#include <meshBase.H>

using namespace nemAux;

class SizeFieldBase
{
  // constructors and destructors
  public:
    SizeFieldBase():mesh(NULL),dev_mult(1.5),maxIsmin(1),sizeFactor(1.)
    {
      std::cout << __FILE__ << __LINE__ << std::endl;
      std::cout << "Size Factor = " << sizeFactor << std::endl;
      std::cout << "SizeFieldBase constructed" << std::endl;
    }
    
    virtual ~SizeFieldBase()
    {
      std::cout << "SizeFieldBase destroyed" << std::endl;
    }

  // SizeFieldBase factory methods
  public:
    // will return derived type GradSizeField|ValSizeField|ErrorSF with SF in mesh's dataSet
    static SizeFieldBase* Create(meshBase* _mesh, std::string method, int arrayID);
    static SizeFieldBase* Create(meshBase* _mesh, std::string method, int arrayID,
                                double _dev_mult, bool _maxIsmin, double _sizeFactor=1.);
    static std::unique_ptr<SizeFieldBase> 
    CreateUnique(meshBase* _mesh, std::string method, int arrayID, double _dev_mult, 
            bool _maxIsmin, double _sizeFactor=1.)
    {
      std::cout << __FILE__ << __LINE__ << std::endl;
      std::cout << "Size Factor = " << _sizeFactor << std::endl;
      return std::unique_ptr<SizeFieldBase>(SizeFieldBase::Create(_mesh,method,arrayID,_dev_mult,_maxIsmin, _sizeFactor));
    }

  // method
  public:
    // compute size field and insert as cell data into mesh's dataSet
    virtual void computeSizeField(int arrayID) = 0;
    void setSizeFactor(double sf) {sizeFactor = sf;};
  
  protected:
    meshBase* mesh;
    double dev_mult;
    bool maxIsmin;
    vtkSmartPointer<vtkDataArray> da;
    std::string sfname;
    double sizeFactor;


  // helpers
  protected:
    // identifies cells to refine and mutates current size values
    // into a compatible size field for the mesh
    void mutateValues(std::vector<double>& values);
    void initialize(meshBase* _mesh, int arrayID, double _dev_mult, 
                    bool _maxIsmin, const std::string& arrName);
};




//class ErrorSF: public SizeFieldBase
//{
//
//}

#endif
