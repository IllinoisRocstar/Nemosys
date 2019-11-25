SetFactory("OpenCASCADE");
Rectangle(1) = {0,0,.1,1,1};
Mesh.CharacteristicLengthMin = 1000;
Mesh.CharacteristicLengthMax = 10000;
Mesh 2;
Save "plane_source_simple.vtk";
