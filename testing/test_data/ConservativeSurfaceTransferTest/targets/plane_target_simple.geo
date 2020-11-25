SetFactory("OpenCASCADE");
Rectangle(1) = {0,0,0,1,1};
Mesh.CharacteristicLengthMin = 0.5;
Mesh.CharacteristicLengthMax = 1;
Mesh 2;
Save "plane_target_simple.vtk";
