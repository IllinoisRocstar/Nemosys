#ifndef ORDEROFACCURACY_H
#define ORDEROFACCURACY_H

#include <meshBase.H>

// f1,f2,f3 -> h1 < h2 < h3

class OrderOfAccuracy
{

  public: 
    OrderOfAccuracy(meshBase* _f1, meshBase* _f2, meshBase* _f3,
                    const std::vector<int>& _arrayIDs);

    std::vector<std::vector<double>> computeOrderOfAccuracy();
    std::vector<std::vector<double>> computeGridConvergenceIndex();

 
  private:
    meshBase *f1, *f2, *f3;
    const std::vector<int> arrayIDs;
    std::vector<int> diffIDs;
    std::vector<int> relEIDs;
    std::vector<std::string> f3ArrNames, f2ArrNames;    
    double r21, r32;
    std::vector<std::vector<double>> diffF3F2;
    std::vector<std::vector<double>> diffF2F1;  
  
    std::vector<std::vector<double>> 
    computeDiffAndRelativeError(meshBase* mesh, const std::vector<std::string>& newArrNames);
};


#endif
