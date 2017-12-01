#ifndef SIZEFIELDGEN_H
#define SIZEFIELDGEN_H

#include <meshBase.H>

class SizeFieldBase
{
  // constructors and destructors
  public:
    SizeFieldBase():mesh(NULL),dev_mult(1.5),maxIsmin(1)
    {
      std::cout << "SizeFieldBase constructed" << std::endl;
    }
    
    virtual ~SizeFieldBase()
    {
      std::cout << "SizeFieldBase destroyed" << std::endl;
    }

  // 'virtual' constructors
  public:
    // will return derived type GradSizeField|ValSizeField|ErrorSF with SF in mesh's dataSet
    static SizeFieldBase* Create(meshBase* _mesh, std::string method, int arrayID);
    static SizeFieldBase* Create(meshBase* _mesh, std::string method, int arrayID,
                                double _dev_mult, bool _maxIsmin);

  // method
  public:
    // compute size field and insert as cell data into mesh's dataSet
    virtual void computeSizeField(int arrayID) = 0;
  
  protected:
    meshBase* mesh;
    double dev_mult;
    bool maxIsmin;

  // helpers
  protected:
    // identifies cells to refine and mutates current size values
    // into a compatible size field for the mesh
    void mutateValues(std::vector<double>& values);
};




//class ErrorSF: public SizeFieldBase
//{
//
//}

#endif
