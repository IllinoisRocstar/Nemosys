#ifndef MESHUSER_H
#define MESHUSER_H
#include <meshBase.H>

class meshUser
{
  public:
    
    meshUser():mesh(){};

    meshUser(std::string _fname): mesh(),fname(_fname)
    {
      mesh = meshBase::Create(fname);
      std::string ext = find_ext(fname);
      // setting write extension for proper vtk file type output
      write_ext = ext; 
      if (ext == ".vtk" || ext == ".vol" || ext == ".msh") // default is vtu
      {
        std::string new_ext(".vtu");
        write_ext = new_ext;
      }  
      
      std::cout << "user constructed" << std::endl;
    }


    
    ~meshUser()
    {
      if (mesh)
      {
        delete mesh;
        mesh = 0;
      }
      std::cout << "user destroyed" << std::endl;
    }

    // information
  public:
    // get number of points in mesh
    int getNumberOfPoints(); 
    // get number of cells in mesh
    int getNumberOfCells();
    // generate a report of the mesh
    void report();      

  // access
  public:
    // get point with id
    std::vector<double> getPoint(int id);
    // get cell with id : returns point indices and respective coordinates
    std::map<int,std::vector<double>> getCell(int id); 
    // get mesh object
    meshBase* getMesh(){ return mesh; };

  // mesh generation
  public:
    int generateMesh(std::string filename, std::string meshEngine); //meshEngine=netgen||gmsh
  
  // data transfer
  public:
    int transfer(meshUser* target, std::string method, int arrayID);
    int transfer(meshUser* target, std::string method);

  // size field generation
  public:
    void generateSizeField(std::string method, int arrayID, double dev_mult, bool maxIsmin);

  // refinement
  public:
    void refineMesh(std::string method, int arrayID, double dev_mult, bool maxIsmin); 

  // diagnostics
  public:
    // print coordinates of point with id
    void printPoint(int id);
    // print point ids and coordinates of cell with id
    void printCell(int id);

  // write and conversion
  public:
    void write(); 
    void write(std::string fname);      
    // convert to gmsh format without data (only volume elements)
    void writeMSH(std::string fname);
    // convert to gmsh format with specified point or cell data
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID); 
  
  private:
    meshBase* mesh;
    std::string fname;
    std::string write_ext;
};

#endif
