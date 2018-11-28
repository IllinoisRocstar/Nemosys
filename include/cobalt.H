#ifndef cobalt_H
#define cobalt_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <vtkCellTypes.h>

// stl
#include <memory>


class meshBase;

namespace COBALT
{

class cobalt
{

    public:
      cobalt(const std::shared_ptr<meshBase> fullMesh, const std::string inFnameVtk, const std::string outFnameCgr, const std::string outFnameCgi);

    public:
      void write() const;

    private:
      std::shared_ptr<meshBase> volMeshBase;
      std::shared_ptr<meshBase> surfMeshBase;
      std::string inFnameVtk;
      std::string outFnameCgr;
      std::string outFnameCgi;
};

}
#endif
