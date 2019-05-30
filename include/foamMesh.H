#ifndef foam_H
#define foam_H

#ifdef HAVE_CFMSH 
#include "meshBase.H"

// vtk
#include "vtkUnstructuredGrid.h"

// OpenFOAM
#include "argList.H"
#include "fvOptions.H"

class meshBase;

namespace FOAM
{

class foamMesh : public meshBase
{

    public:
      foamMesh(bool readDB=false);
      ~foamMesh();

  // --- diagnostics
  public:
    // generate a report of the mesh
    virtual void report(); 

    // access
    public:
      // general purpose mesh information read method
      void read(const std::string& fname);
      // get point with id
      std::vector<double> getPoint(int id) const;
      // get 3 vecs with x,y and z coords
      std::vector<std::vector<double>> getVertCrds() const;
      // get cell with id : returns point indices and respective coordinates
      std::map<int, std::vector<double>> getCell(int id) const;
      std::vector<std::vector<double>> getCellVec(int id) const;
      // get diameter of circumsphere of each cell
      std::vector<double> getCellLengths() const;
      // get center of a cell
      std::vector<double> getCellCenter(int cellID) const;
      std::vector<double> getCellCenter(int cellID);
      // get cell type as an integer
      // assumes all elements are the same type
      int getCellType() const;
      // get edge lengths of dataSet
      void inspectEdges(const std::string &ofname) const;
      std::vector<int> getConnectivities() const;
      // write the mesh to file named fname
      void write(const std::string &fname) const;


    // processing
    public:
       vtkSmartPointer<vtkDataSet> extractSurface();
    
    // internal methods
    private:
      void genMshDB();
      void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
              const int cellType, std::vector<int>& vrtIds);


    // openfoam data structure 
    private:
        Foam::argList* _args;
        Foam::Time* _runTime;
        Foam::fvMesh* _fmesh;
    
    // internal data
    private:
      std::shared_ptr<meshBase> _volMB;
};

}

#endif

#endif
