#include "diffMesh.H"

#include <cmath>

/*
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
*/

namespace NEM {
namespace MSH {

bool is_close(double val1, double val2, double floor, double rel_tol) {
  if (std::fabs(val1) <= floor && std::fabs(val2) <= floor) {
    return true;
  } else {
    return std::fabs((val1 - val2) / val1) < rel_tol;
  }
}

int diffMesh(geoMeshBase *gmb1, geoMeshBase *gmb2, double floor, double relTol,
             double numCellsTol, double numPointsTol) {
  // Check if equal number of points and cells.
  {
    auto points1 = gmb1->getNumberOfPoints();
    auto points2 = gmb2->getNumberOfPoints();
    if (numPointsTol > 0) {
      if (!is_close(points1, points2, 0, numPointsTol)) {
        std::cerr << "Meshes don't have similar number of points: " << points1
                  << " vs " << points2 << std::endl;
        return 1;
      }
    } else {
      if (points1 != points2) {
        std::cerr << "Meshes don't have same number of points: " << points1
                  << " vs " << points2 << std::endl;
        return 1;
      }
    }
  }
  {
    auto cells1 = gmb1->getNumberOfCells();
    auto cells2 = gmb2->getNumberOfCells();
    if (numCellsTol > 0) {
      if (!is_close(cells1, cells2, 0, numCellsTol)) {
        std::cerr << "Meshes don't have similar number of cells: " << cells1
                  << " vs " << cells2 << std::endl;
        return 1;
      }
    } else {
      if (cells1 != cells2) {
        std::cerr << "Meshes don't have same number of cells: " << cells1
                  << " vs " << cells2 << std::endl;
        return 1;
      }
    }
  }

  if (numPointsTol <= 0) {
    // Check if point coordinates match within tolerance
    for (int i = 0; i < gmb1->getNumberOfPoints(); ++i) {
      std::array<double, 3> coord1{};
      std::array<double, 3> coord2{};
      gmb1->getPoint(i, &coord1);
      gmb2->getPoint(i, &coord2);
      for (int j = 0; j < 3; ++j) {
        if (!is_close(coord1[j], coord2[j], floor, relTol)) {
          std::cerr << "Meshes differ in point coordinates" << std::endl;
          std::cerr << "Index " << i << " Component " << j << std::endl;
          std::cerr << "Coord 1 " << std::setprecision(15) << coord1[j]
                    << " Coord 2 " << std::setprecision(15) << coord2[j]
                    << std::endl;
          std::cerr << "Meshes differ in point coordinates" << std::endl;
          return 1;
        }
      }
    }
  }

  if (numCellsTol <= 0) {
    for (int i = 0; i < gmb1->getNumberOfCells(); ++i) {
      vtkSmartPointer<vtkPoints> points1 = gmb1->getCell(i)->GetPoints();
      vtkSmartPointer<vtkPoints> points2 = gmb2->getCell(i)->GetPoints();

      if (points1->GetNumberOfPoints() != points2->GetNumberOfPoints()) {
        std::cerr << "Meshes differ in cells" << std::endl;
        return 1;
      }

      for (int j = 0; j < points1->GetNumberOfPoints(); ++j) {
        std::vector<double> coord1(3);
        std::vector<double> coord2(3);
        points1->GetPoint(j, coord1.data());
        points2->GetPoint(j, coord2.data());
        for (int k = 0; k < 3; ++k) {
          if (!is_close(coord1[k], coord2[k], floor, relTol)) {
            std::cerr << "Meshes differ in cells" << std::endl;
            return 1;
          }
        }
      }
    }
  }

  int numPointArr1 = gmb1->getNumberOfPointDataArrays();
  int numPointArr2 = gmb2->getNumberOfPointDataArrays();

  if (numPointArr1 != numPointArr2) {
    std::cerr << "Meshes have different numbers of point data" << std::endl;
    std::cerr << "Mesh 1 has " << numPointArr1 << std::endl;
    std::cerr << "Mesh 2 has " << numPointArr2 << std::endl;
    return 1;
  }

  for (int i = 0; i < numPointArr1; ++i) {
    auto arr1 = gmb1->getPointDataArrayCopy(i);
    auto arr2 = gmb2->getPointDataArrayCopy(arr1->GetName());
    if (!arr2) {
      std::cerr << "Mesh 2 does not have a point data array of name "
                << arr1->GetName() << std::endl;
      return 1;
    }
    /*
    if (arr1->GetDataType() != arr2->GetDataType()) {
      std::cerr << "For point data array " << arr1->GetName()
                << ",\narrays do not have same data type:\n"
                << arr1->GetDataTypeAsString() << " vs "
                << arr2->GetDataTypeAsString() << std::endl;
      return 1;
    }
    */
    if (arr1->GetNumberOfComponents() != arr2->GetNumberOfComponents()) {
      std::cerr << "For point data array " << arr1->GetName()
                << ",\narrays do not have same number of components:\n"
                << arr1->GetNumberOfComponents() << " vs "
                << arr2->GetNumberOfComponents() << std::endl;
      return 1;
    }
    if (numPointsTol <= 0) {
      auto da1 = vtkDataArray::FastDownCast(arr1);
      auto da2 = vtkDataArray::FastDownCast(arr2);
      if (da1 && da2) {
        int numComponents = da1->GetNumberOfComponents();
        for (int j = 0; j < gmb1->getNumberOfPoints(); ++j) {
          std::vector<double> comps1(numComponents);
          std::vector<double> comps2(numComponents);
          da1->GetTuple(j, comps1.data());
          da2->GetTuple(j, comps2.data());
          for (int k = 0; k < numComponents; ++k) {
            if (!is_close(comps1[k], comps2[k], floor, relTol)) {
              std::cerr << "For point data array " << da1->GetName()
                        << std::endl;
              std::cerr << "Meshes differ in point data values at point " << j
                        << " component " << k << std::endl;
              std::cerr << comps1[k] << " " << comps2[k] << std::endl;
              return 1;
            }
          }
        }
      }
    }
  }

  int numCellArr1 = gmb1->getNumberOfCellDataArrays();
  int numCellArr2 = gmb2->getNumberOfCellDataArrays();

  if (numCellArr1 != numCellArr2) {
    std::cerr << "Meshes have different numbers of cell data" << std::endl;
    std::cerr << "Mesh 1 has " << numCellArr1 << std::endl;
    std::cerr << "Mesh 2 has " << numCellArr2 << std::endl;
    return 1;
  }

  for (int i = 0; i < numCellArr1; ++i) {
    auto arr1 = gmb1->getCellDataArrayCopy(i);
    if (arr1->GetName() == gmb1->getGeoEntArrayName()) {
      continue;
    }
    auto arr2 = gmb2->getCellDataArrayCopy(arr1->GetName());
    if (!arr2) {
      std::cerr << "Mesh 2 does not have a cell data array of name "
                << arr1->GetName() << std::endl;
      return 1;
    }
    /*
    if (arr1->GetDataType() != arr2->GetDataType()) {
      std::cerr << "For cell data array " << arr1->GetName()
                << ",\narrays do not have same data type:\n"
                << arr1->GetDataTypeAsString() << " vs "
                << arr2->GetDataTypeAsString() << std::endl;
      return 1;
    }
    */
    if (arr1->GetNumberOfComponents() != arr2->GetNumberOfComponents()) {
      std::cerr << "For cell data array " << arr1->GetName()
                << ",\narrays do not have same number of components:\n"
                << arr1->GetNumberOfComponents() << " vs "
                << arr2->GetNumberOfComponents() << std::endl;
      return 1;
    }
    if (numCellsTol <= 0) {
      auto da1 = vtkDataArray::FastDownCast(arr1);
      auto da2 = vtkDataArray::FastDownCast(arr2);
      if (da1 && da2) {
        int numComponents = da1->GetNumberOfComponents();
        for (int j = 0; j < gmb1->getNumberOfCells(); ++j) {
          std::vector<double> comps1(numComponents);
          std::vector<double> comps2(numComponents);
          da1->GetTuple(j, comps1.data());
          da2->GetTuple(j, comps2.data());
          for (int k = 0; k < numComponents; ++k) {
            if (!is_close(comps1[k], comps2[k], floor, relTol)) {
              std::cerr << "For cell data array " << da1->GetName()
                        << std::endl;
              std::cerr << "Meshes differ in cell data values at cell " << j
                        << " component " << k << std::endl;
              std::cerr << comps1[k] << " " << comps2[k] << std::endl;
              return 1;
            }
          }
        }
      }
    }
  }
  std::cout << "Meshes are the same" << std::endl;
  return 0;
}

}  // namespace MSH
}  // namespace NEM
