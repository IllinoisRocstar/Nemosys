//Plane Definition
SetFactory("OpenCASCADE");
z = 0.800000;
pl = 2.000000;
lc = 0.01;
Point(1) = {pl, 0, z, lc};
Point(2) = {0, 0, z, lc};
Point(3) = {0, pl, z, lc};
Point(4) = {pl,  pl, z, lc} ;
Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;
Line Loop(1) = {1,2,3,4} ;
Plane Surface(1) = {1} ;
Periodic Line {1} = {-3} ;
Periodic Line {2} = {-4} ;
//All Spheres
lc = 0.01;
SetFactory("OpenCASCADE");
Sphere(5) = {1.723565,0.958345,0.437759,0.568406};
Sphere(6) = {-0.276435,0.958345,0.437759,0.568406};
Sphere(7) = {1.187667,1.958345,0.365962,0.568406};
Sphere(8) = {1.187667,-0.041655,0.365962,0.568406};
Sphere(9) = {1.651768,0.422446,1.437759,0.568406};
Sphere(10) = {-0.348232,0.422446,1.437759,0.568406};
Sphere(11) = {-0.348232,0.422446,-0.562241,0.568406};
Sphere(12) = {1.651768,0.422446,-0.562241,0.568406};
Sphere(13) = {0.723565,0.886548,1.901861,0.568406};
Sphere(14) = {0.723565,0.886548,-0.098139,0.568406};
Sphere(15) = {0.651768,1.886548,1.365962,0.568406};
Sphere(16) = {0.651768,-0.113452,1.365962,0.568406};
Physical Volume("sphere_1.00") = {5:16};
s0 = BooleanIntersection{ Surface{1}; }{Volume{5:16}; };
Physical Surface("sphere_1.00") = s0();
d = BooleanDifference{ Surface{1}; Delete;}{Volume{5:16}; Delete;};
Physical Surface("Domain") = d();
//Meshhhh
Mesh.Algorithm = 6; //Frontal
Mesh.Optimize = 1; // Gmsh smoother, works with boundary layers (netgen version does not).
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthMax = 0.05;
Mesh.RemeshAlgorithm = 1; //(0) nosplit (1) automatic (2) split only with metis
Mesh 2;
Save StrCat(StrPrefix(General.FileName), ".msh");
