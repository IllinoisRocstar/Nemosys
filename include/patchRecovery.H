#ifndef PATCHRECOVERY_H
#define PATCHRECOVERY_H

#include <orthoPoly3D.H>
#include <Cubature.H>

class PatchRecovery
{
  
  public:
    PatchRecovery(meshBase* nodeMesh, int _order, const std::vector<int>& arrayIDs);
    ~PatchRecovery(){} 
    void extractAxesAndData(pntDataPairVec& pntsAndData, 
														std::vector<std::vector<double>>& coords,
                            std::vector<VectorXd>& data,
                            const std::vector<int>& numComponents,
                            int& pntNum);

    void recoverNodalSolution();

  private:
    std::unique_ptr<GaussCubature> cubature;     
    std::vector<orthoPoly3D> polyApproximants;
		int order;



    // disabled
    PatchRecovery(const PatchRecovery& that) = delete;
    PatchRecovery& operator=(const PatchRecovery& that) = delete; 
};






#endif
