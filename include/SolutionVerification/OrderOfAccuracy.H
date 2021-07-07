#ifndef NEMOSYS_ORDEROFACCURACY_H_
#define NEMOSYS_ORDEROFACCURACY_H_

#include "nemosys_export.h"

#include "Mesh/meshBase.H"

#include <vector>

/**
 * The OrderOfAccuracy class computes the grid convergence
 * index (GCI) given three successively refined meshes. The
 * GCI is particularly robust since it accounts for the order
 * of convergence as well as the grid refinement ratio.
 *
 * If a vector field is being analyzed, a GCI will be generated
 * for each component.
 *
 * For
 * more information see:
 * Roache, P.J.
 * Verification and Validation in Compuational Science and
 * Engineering, Hermosa Publishers, Albuquerque, New Mexico,
 * 1998
 */
class NEMOSYS_EXPORT OrderOfAccuracy {
 public:
  OrderOfAccuracy(meshBase *_f3, meshBase *_f2, meshBase *_f1,
                  std::vector<int> _arrayIDs,
                  std::string transferType = "Consistent Interpolation",
                  double targetGCI = 1.1);

  ~OrderOfAccuracy() = default;

  /**
   * @brief Compute order of accuracy.
   * @return Order of accuracy for each component of each field
   * based on errors between the coarse and finer mesh and the
   * finer and finest mesh.
   */
  std::vector<std::vector<double>> computeOrderOfAccuracy();

  /**
   * @brief Compute the GCI with respect to the finer and finest
   * grids.
   * @return GCI for each component of each field transferred
   * from the finer mesh (f2) to the finest mesh (f1).
   */
  std::vector<std::vector<double>> computeGCI_21();

  /**
   * @brief Compute the GCI with respect to the coarsest and
   * finer grids.
   * @return GCI for each component of each field transferred
   * from the coarse mesh (f3) to the finer mesh (f2).
   */
  std::vector<std::vector<double>> computeGCI_32();

  /**
   * @brief Compute the grid resolution required to obtain the
   * desired gci/level of accuracy for each component across
   * all mesh fields. Given a vector field, select the maximum
   * refinement ratio in the components to ensure the desired
   * GCI is met across all components.
   * @return Necessary grid refinement ratio for each component
   * of each field.
   */
  std::vector<std::vector<double>> computeResolution(double gciStar);
  void computeRichardsonExtrapolation();

  /**
   * @brief Returns order of accuracy for each component of each field
   * @return Order of accuracy of each component of each field
   */
  std::vector<std::vector<double>> getOrderOfAccuracy() const {
    return orderOfAccuracy;
  }

  /**
   * @brief Check if grids are in asymptotic range of convergence based
   * on the computed GCIs (GCI_32 and GCI_21). We are within asymptotic
   * range if all values are close to 1
   * @return Weighted ratio of GCIs for each component of each field in grids
   */
  std::vector<std::vector<double>> checkAsymptoticRange();

  /**
   * @brief Refine the coarsest mesh based on the target GCI, write
   * to file with name ofname.
   * @param gciStar target GCI
   * @param ofname refined mesh filename
   */
  void computeMeshWithResolution(double gciStar, const std::string &ofname);

  std::vector<std::vector<double>> computeDiff(
      meshBase *mesh, const std::vector<std::string> &newArrNames);

  std::vector<std::vector<double>> getDiffF2F1() const { return diffF2F1; }

  double getTargetGCI() const { return targetGCI; }

 private:
  meshBase *f3, *f2, *f1;  ///< meshes from least to most refined (f3 -
                           ///< coarsest, f2 - finer, f1 - finest)
  const std::vector<int>
      arrayIDs;  ///< array ids for arrays to be considered in analysis
  std::vector<int>
      diffIDs;  ///< array ids for differences^2 in solutions between meshes
  std::vector<int>
      relEIDs;  ///< array ids for integral of solutions on the most refined
                ///< mesh, used in computing relative discretization error
  std::vector<int> realDiffIDs;  ///< array ids for actual difference data
  std::vector<std::string> f3ArrNames,
      f2ArrNames;   ///< names of data arrays transferred from f3 to f2 and from
                    ///< f2 to f1
  double r21, r32;  ///< effective grid refinement ratios
  std::vector<std::vector<double>>
      diffF3F2;  ///< L2 norm of the difference in solution between f3 and f2
  std::vector<std::vector<double>>
      diffF2F1;  ///< L2 norm of the difference in solution between f2 and f1
  std::vector<std::vector<double>>
      GCI_32;  ///< GCIs with respect to the coarse mesh (f3)
  std::vector<std::vector<double>>
      GCI_21;  ///< GCIs with respect to the finer mesh (f2)
  std::vector<std::vector<double>>
      orderOfAccuracy;  ///< Order of accuracy (p) with respect to refinements
                        ///< f3-f2-f1
  double targetGCI;
};

#endif  // NEMOSYS_ORDEROFACCURACY_H_
