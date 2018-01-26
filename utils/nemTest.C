#include <meshBase.H>

int main(int argc, char* argv[])
{
  
  meshBase* mesh = meshBase::Create("refined1_transfer.vtu");
  mesh->getIntegrationPointsAtCell(50);
  delete mesh;
  return 0;

}
