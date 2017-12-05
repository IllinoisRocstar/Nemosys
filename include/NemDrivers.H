#ifndef NEMDRIVERS_H
#define NEMDRIVERS_H

#include <meshUser.H>


/* ---------------------------------------------------- 
TODO: This class will contain drivers for several utils:
    1) Mesh Generation :
      - ng_meshing_params
      - stl file
    2) Data Transfer:
      - mesh file with data (source)
      - mesh file without data (target)
      - name of array to transfer
        OR no name (transfers all)
      - method to use for transfer (currently only "FE")
    3) Adaptive Mesh Refinement:
      - mesh file with data
      - method of refinement ("gradient", "value", "error")
      - name of array to use for size field gen (w/ method)
    4) Uniform Mesh Refinement:
      - mesh file w/ or w/o data
      - scaling factor for edges   
    5) Format conversion:
      - vt* to .msh, .vol etc.
      - .msh to vtu
      - .vol to vtu 

   -------------------------------------------------- */



class TransferDriver
{
  public:
    TransferDriver(std::string srcmsh, std::string trgmsh,
                   std::string ofname)
    {
      source = new meshUser(srcmsh);
      target = new meshUser(trgmsh);
      Timer T;
      T.start();
      source->transfer(target, "FE");
      T.stop();
      std::cout << "Time spent transferring data (ms)" << T.elapsed() << std::endl;
      target->write(ofname); 
    }

    TransferDriver(std::string srcmsh, std::string trgmsh, int arrayID,
                   std::string ofname)
    {
      source = new meshUser(srcmsh);
      target = new meshUser(trgmsh);
      Timer T;
      T.start();
      source->transfer(target, "FE", arrayID);
      T.stop();
      std::cout << "Time spent transferring data (ms)" << T.elapsed() << std::endl;
      target->write(ofname); 
    }

    TransferDriver(std::string srcmsh, std::string trgmsh, std::vector<int> arrayIDs,
                   std::string ofname)
    {
      source = new meshUser(srcmsh);
      target = new meshUser(trgmsh);
      Timer T;
      T.start();
      source->transfer(target, "FE", arrayIDs);
      T.stop();
      std::cout << "Time spent transferring data (ms)" << T.elapsed() << std::endl;
      target->write(ofname); 
    }

    ~TransferDriver()
    {
      if (source)
      {
        delete source;
        source = 0;
      }
      if (target)
      {
        delete target;
        target = 0;
      }
    }

  private:
    meshUser* source;
    meshUser* target; 

};

#endif
