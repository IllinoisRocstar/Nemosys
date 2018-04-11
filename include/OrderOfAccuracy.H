#ifndef ORDEROFACCURACY_H
#define ORDEROFACCURACY_H

#include <meshBase.H>

// f1,f2,f3 -> h1 < h2 < h3

class OrderOfAccuracy
{

  public: 
    OrderOfAccuracy(meshBase* _f1, meshBase* _f2, meshBase* _f3,
                    const std::vector<int>& _arrayIDs);

    ~OrderOfAccuracy();

    // get order of accuracy by solving for p in (f3-f2)/(r23^p - 1) = (r21^p*(f2-f1)/(r21^p-1)
    std::vector<std::vector<double>> computeOrderOfAccuracy();
    // get the grid convergence index for the mesh f2
    std::vector<std::vector<double>> computeGCI_21();
    // get the grid convergence index for the mesh f3
    std::vector<std::vector<double>> computeGCI_32();
    // get the grid resolution required to obtain the desired gci/level of accuracy
    std::vector<std::vector<double>> computeResolution(double gciStar);
    // return order of accuracy
    std::vector<std::vector<double>> getOrderOfAccuracy() { return orderOfAccuracy; };
    // check if the grids are in the asymptotic range of convergnece
    // they must be in order for the analysis performed by this class to make sense
    std::vector<std::vector<double>> checkAsymptoticRange();
    // compute the difference in solutions specified by arrayIDs
    // the arrays have been transfered to 'mesh' and the L2 norm of the difference is returned
    // if realDiff is on, the actual difference (not sqr) is added to 'mesh' as pointdata 
    std::vector<std::vector<double>> 
    computeDiff(meshBase* mesh, const std::vector<std::string>& newArrNames);
    void computeRichardsonExtrapolation();
    void computeMeshWithResolution(double gciStar, const std::string& ofname);
 
  private:
    // meshes from most to least refined
    meshBase *f1, *f2, *f3;
    // arrays to consider for analysis
    const std::vector<int> arrayIDs;
    // pointdata array ids for differences in solutions between meshes
    std::vector<int> diffIDs;
    // pointdata array ids for integral of the solutions on the most refined mesh
    // this is used to compute the relative discretization error
    std::vector<int> relEIDs;
    // array ids for actual difference data
    std::vector<int> realDiffIDs;
    // names of data arrays transfered from f3 to f2 and from f2 to f1
    std::vector<std::string> f3ArrNames, f2ArrNames;    
    // effective grid refinement ratios
    double r21, r32;
    // L2 norm of the difference in solution between f3 and f2
    std::vector<std::vector<double>> diffF3F2;
    // L2 norm of the difference in solution between f2 and f1
    std::vector<std::vector<double>> diffF2F1; 
    // grid convergence index for the coarse mesh 
    std::vector<std::vector<double>> GCI_32;
    // grid convergence index for the middle refined mesh
    std::vector<std::vector<double>> GCI_21;
    // observed order of accuracy for solution specified by arrayIDs
    std::vector<std::vector<double>> orderOfAccuracy; 

};


#endif
