#ifndef SIZEFIELDGEN_H
#define SIZEFIELDGEN_H

#include <meshBase.H>

class SizeFieldGen
{
  // constructors and destructors
  public:
    SizeFieldGen():mesh(NULL),dev_mult(1.5),maxIsmin(1)
    {
      std::cout << "SizeFieldGen constructed" << std::endl;
    }
    
    virtual ~SizeFieldGen()
    {
      std::cout << "SizeFieldGen destroyed" << std::endl;
    }

  // 'virtual' constructors
  public:
    // will return derived type GradSF|ValSF|ErrorSF with SF in mesh's dataSet
    static SizeFieldGen* Create(meshBase* _mesh, std::string method, int arrayID);
    static SizeFieldGen* Create(meshBase* _mesh, std::string method, int arrayID,
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

class GradSF: public SizeFieldGen
{
  // constructors and destructors
  public:
    GradSF(meshBase* _mesh, int arrayID, double _dev_mult, bool _maxIsmin);
    ~GradSF()
    {
      std::cout << "GradSF destroyed" << std::endl;
    }
    
  // methods
  public:
    // computes the gradient of point data at a cell using 
    // derivatives of shape interpolation functions
    std::vector<double> computeGradAtCell(int cellID, int arrayID);
    // compute L2 norm of gradient of point data at each cell 
    std::vector<double> computeL2GradAtAllCells(int arrayID);
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID){}

  private:
    vtkDataArray* da;
    std::string sfname; 
};

class ValSF: public SizeFieldGen
{

  // constructors and destructors
  public:
    ValSF(meshBase* _mesh, int arrayID, double _dev_mult, bool _maxIsmin);
    ~ValSF()
    {
      std::cout << "ValSF destroyed" << std::endl;
    }

  // methods
  public:
    // computes value of point data at a cell center using shape interpolation functions
    std::vector<double> computeValAtCell(int cellID, int arrayID);
    // compute value of point data at center of each cell
    std::vector<std::vector<double>> computeValAtAllCells(int arrayID);
    // compute L2 norm of value of point data at center of each cell
    std::vector<double> computeL2ValAtAllCells(int array);
    // compute size field and insert as cell data into mesh's dataSet
    void computeSizeField(int arrayID);

  private:
    vtkDataArray* da;
    std::string sfname;
};

//class ErrorSF: public SizeFieldGen
//{
//
//}

#endif
