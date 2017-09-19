#include <MAdLib.h>

int main(int argc, char* argv[])
{
  MAd::pGModel gmodel = NULL;
  MAd::GM_create(&gmodel, "GeoModel");
  MAd::GM_read(gmodel, argv[1]);
  
  MAd::pMesh mesh = MAd::M_new(gmodel);
  MAd::M_load(mesh, argv[2]); 
  
  MAd::PWLSField* sizeField = new MAd::PWLSField(mesh);
  
  sizeField->setCurrentSize();
  sizeField->scale(0.5);
  
  MAd::MeshAdapter* adapter = new MAd::MeshAdapter(mesh, sizeField);
  adapter->run();

  MAd::M_writeMsh(mesh, "newMesh.msh", 2);

  return 0;
}
