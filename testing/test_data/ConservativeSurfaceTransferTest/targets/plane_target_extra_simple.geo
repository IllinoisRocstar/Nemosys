SetFactory("OpenCASCADE");
Rectangle(1) = {0,0,0,1,1};
Mesh.CharacteristicLengthMin = 5;
Mesh.CharacteristicLengthMax = 10;
Mesh 2;
Save "plane_target_extra_simple.vtk";
