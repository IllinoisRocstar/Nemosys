#ifndef SURFACELAMBDAMUSMOOTH_H
#define SURFACELAMBDAMUSMOOTH_H

#include "MeshedSurfaces.H"
#include "Time.H"
#include "argList.H"
#include "boundBox.H"
#include "edgeMesh.H"
#include "fileName.H"
#include "fvMesh.H"
#include "matchPoints.H"

class surfaceLambdaMuSmooth {
 public:
  /** @brief Constructor */
  surfaceLambdaMuSmooth();

  /** @brief Destructor */
  ~surfaceLambdaMuSmooth();

  /**
   * @brief smooths out surface features using tolerances
   * @param surfFileName Input surface file name
   * @param outFileName Output surface file name
   * @param lambda Lambda tolerance value
   * @param mu     Mu tolerance value
   * @param iters  Number of iterations allowed
   * @param addFtrFl To include feature file or not
   */
  void execute(const Foam::fileName surfFileName,
               const Foam::fileName outFileName, const Foam::scalar lambda,
               const Foam::scalar mu, const Foam::label iters,
               const bool addFtrFl);

 private:
  Foam::tmp<Foam::pointField> avg(const Foam::meshedSurface &s,
                                  const Foam::bitSet &fixedPoints);

  void getFixedPoints(const Foam::edgeMesh &feMesh,
                      const Foam::pointField &points,
                      Foam::bitSet &fixedPoints);
};

#endif  // SURFACELAMBDAMUSMOOTH_H
