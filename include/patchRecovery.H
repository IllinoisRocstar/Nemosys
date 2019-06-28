#ifndef PATCHRECOVERY_H
#define PATCHRECOVERY_H

#include "nemosys_export.h"

#include "Cubature.H"


// pair type for coordinate and data there
//typedef std::pair<std::vector<double>, std::vector<double>> pntDataPair;

// holds gauss points and data at these points as pairs
//typedef std::vector<pntDataPair> pntDataPairVec;

class NEMOSYS_EXPORT PatchRecovery
{
  public:
    PatchRecovery(meshBase *nodeMesh, int _order,
                  const std::vector<int> &arrayIDs);

    ~PatchRecovery() = default;

    // disabled copy and copy-assignment constructors
    PatchRecovery(const PatchRecovery &that) = delete;
    PatchRecovery &operator=(const PatchRecovery &that) = delete;

    // approximates nodal values with polynomial approximant over patch of node
    // sets recovered nodal values as new point data
    void recoverNodalSolution(bool ortho);
    /* 1) approximates nodal values with polynomial approximant over patch of node
       2) sets error^2 in recovered values as new point data
       3) sets RMSE integrated over each cell, normalized by cell vol as new cell data
       4) sets patch averaged element size as new point data
       returns array of integrals of errors over entire mesh */
    std::vector<std::vector<double>> computeNodalError();

  private:
    std::unique_ptr<GaussCubature> cubature;
    int order;

    // get min max of each coord (for regularization)
    std::vector<double>
    getMinMaxCoords(const std::vector<std::vector<double>> &coords) const;
    // regularize coordinates as x_new = -1 + 2*(x_old - x_min)/(x_max - x_min)
    // for proper conditioning of least squares matrix
    void regularizeCoords(std::vector<std::vector<double>> &coords,
                          std::vector<double> &genNodeCoord) const;

    // extract coordinates and data from pntDataPair
    void extractAxesAndData(const pntDataPairVec &pntsAndData,
                            std::vector<std::vector<double>> &coords,
                            std::vector<VectorXd> &data,
                            const std::vector<int> &numComponents,
                            int &pntNum) const;
};

#endif
