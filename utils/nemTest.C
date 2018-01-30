#include <Cubature.H>
int main(int argc, char* argv[])
{
  
  meshBase* mesh = meshBase::Create("refined1_transfer.vtu");
  //mesh->getIntegrationPointsAtCell(50);
  GaussCubature* cuby = new  GaussCubature(mesh);
  std::vector<int> tmp = {0,2,3,7};
  cuby->interpolateToGaussPoints(tmp);
	cuby->writeGaussMesh();
  meshBase* mesh1 = meshBase::Create("refined1_transferGaussPoints.vtp");
  mesh1->write("tmptmptmp.vtp");
  delete cuby;
  delete mesh;
  delete mesh1;
  return 0;

}
