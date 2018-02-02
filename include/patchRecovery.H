#ifndef PATCHRECOVERY_H
#define PATCHRECOVERY_H

#include <orthoPoly3D.H>
#include <Cubature.H>

class PatchRecovery
{
  
  public:
    PatchRecovery(meshBase* nodeMesh, const std::vector<int>& arrayIDs);
     

  private:
    std::unique_ptr<GaussCubature> cubature;     
    std::vector<orthoPoly3D> polyApproximants;

    void recoverNodalSolution(const std::vector<int>& arrayIDs);

    // disabled
    PatchRecovery(const PatchRecovery& that) = delete;
    PatchRecovery& operator=(const PatchRecovery& that) = delete; 
};






#endif
