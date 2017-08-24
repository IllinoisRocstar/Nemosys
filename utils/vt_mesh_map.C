#include <iostream>
#include <vtkAnalyzer.H>


int main(int argc, char* argv[])
{
  if (argc < 3)
	{
		std::cerr << " VTK 3D-2D mesh transfer utility" << std::endl
							<< " Usage: " << argv[0] << " VolMeshFile PlaneMeshFile" 
							<< std::endl;
		exit(1);
	}

	// Reading files from command line
	vtkAnalyzer* VolMesh;
	vtkAnalyzer* PlaneMesh;
	VolMesh = new vtkAnalyzer(argv[1]);
	PlaneMesh = new vtkAnalyzer(argv[2]);
	VolMesh->read();
	PlaneMesh->read();
	
	// Report mesh statistics
	VolMesh->report();
	std::cout << std::endl;
	PlaneMesh->report();
	
	/*unsigned long num_points = PlaneMesh->getNumberofPoints();
	
	for (int i = 0; i < num_points; ++i) {
		double* coord = PlaneMesh->getPointCoords(i);
			

		std::cout << coords[0] << " " << coords[1] << " "
							<< coords[2] << std::endl;
	}	*/

	int numComponent;	
	unsigned long num_cells = PlaneMesh->getNumberOfCells();
	for (int i = 0; i < num_cells; ++i) {
	//	PlaneMesh->getCellCoords(i);}
		double** cellCoords = PlaneMesh->getCellCoords(i,numComponent);
		for (int j = 0; j < numComponent ; ++j) {
			std::cout << cellCoords[j][0] << " " << cellCoords[j][1] 
								<< " " << cellCoords[j][2] << std::endl;
			free(cellCoords[j]);
			}
		std::cout << std::endl;
		free(cellCoords);
	}

	delete VolMesh;
	delete PlaneMesh;
	return 0;
}

