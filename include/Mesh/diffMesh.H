#ifndef NEMOSYS_DIFFMESH_H_
#define NEMOSYS_DIFFMESH_H_

#include "nemosys_export.h"
#include "Mesh/geoMeshBase.H"

namespace NEM {
namespace MSH {

/**
 * Compare two @c geoMeshBase objects.
 * @param gmb1 Reference @c geoMeshBase
 * @param gmb2 Other @c geoMeshBase
 * @param floor If a floating point value is less than or equal to @p floor, it
 * is treated as 0. Should be non-negative
 * @param relTol Relative tolerance for comparing floating point values
 * @param numCellsTol Tolerance for number of cells, if positive. If non-positive,
 * exact match required for number of cells and cell data compared at each cell.
 * @param numPointsTol Tolerance for number of points, if positive. If
 * non-positive, exact match required for number of points and point data
 * compared at each cell.
 * @return 0 if @p gmb1 and @p gmb2 are similar; 1 otherwise
 */
NEMOSYS_EXPORT int diffMesh(geoMeshBase *gmb1, geoMeshBase *gmb2,
                            double floor = 1e-9, double relTol = 1e-6,
                            double numCellsTol = -1.,
                            double numPointsTol = -1.);

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_DIFFMESH_H_
