SetFactory("OpenCASCADE");
Rectangle(1) = {0,0,0,1,1};
Mesh.CharacteristicLengthMin = 0.02;
Mesh.CharacteristicLengthMax = 0.05;
Mesh 2;
Save "plane_source.vtk";
