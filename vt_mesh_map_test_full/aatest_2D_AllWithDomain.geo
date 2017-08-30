//Plane Definition
SetFactory("OpenCASCADE");
z = 0.800000;
pl = 1.000000;
lc = 0.01;
Point(1) = {pl, -pl, z, lc};
Point(2) = {-pl, -pl, z, lc};
Point(3) = {-pl, pl, z, lc};
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
Sphere(7) = {-0.276435,-1.041655,0.437759,0.568406};
Sphere(8) = {-0.276435,-1.041655,-1.562241,0.568406};
Sphere(9) = {-0.276435,-1.041655,0.437759,0.568406};
Sphere(10) = {1.723565,-1.041655,0.437759,0.568406};
Sphere(11) = {1.723565,-1.041655,-1.562241,0.568406};
Sphere(12) = {1.723565,0.958345,-1.562241,0.568406};
Sphere(13) = {1.187667,1.958345,0.365962,0.568406};
Sphere(14) = {-0.812333,1.958345,0.365962,0.568406};
Sphere(15) = {-0.812333,-0.041655,0.365962,0.568406};
Sphere(16) = {1.187667,-0.041655,0.365962,0.568406};
Sphere(17) = {1.651768,0.422446,1.437759,0.568406};
Sphere(18) = {-0.348232,0.422446,1.437759,0.568406};
Sphere(19) = {-0.348232,0.422446,-0.562241,0.568406};
Sphere(20) = {1.651768,0.422446,-0.562241,0.568406};
Sphere(21) = {0.723565,0.886548,1.901861,0.568406};
Sphere(22) = {-1.276435,0.886548,1.901861,0.568406};
Sphere(23) = {-1.276435,-1.113452,1.901861,0.568406};
Sphere(24) = {-1.276435,-1.113452,-0.098139,0.568406};
Sphere(25) = {0.723565,-1.113452,1.901861,0.568406};
Sphere(26) = {0.723565,-1.113452,-0.098139,0.568406};
Sphere(27) = {0.723565,0.886548,-0.098139,0.568406};
Sphere(28) = {0.651768,1.886548,1.365962,0.568406};
Sphere(29) = {-1.348232,1.886548,1.365962,0.568406};
Sphere(30) = {-1.348232,-0.113452,1.365962,0.568406};
Sphere(31) = {-1.348232,-0.113452,-0.634038,0.568406};
Sphere(32) = {-1.348232,-0.113452,1.365962,0.568406};
Sphere(33) = {0.651768,-0.113452,1.365962,0.568406};
Sphere(34) = {0.651768,-0.113452,-0.634038,0.568406};
Sphere(35) = {0.651768,1.886548,-0.634038,0.568406};
Physical Volume("sphere_1.00") = {5:35};
s0 = BooleanIntersection{ Surface{1}; }{Volume{5:35}; };
Physical Surface("sphere_1.00") = s0();
d = BooleanDifference{ Surface{1}; Delete;}{Volume{5:35}; Delete;};
Physical Surface("Domain") = d();
//Meshhhh
Mesh.Algorithm = 6; //Frontal
Mesh.Optimize = 1; // Gmsh smoother, works with boundary layers (netgen version does not).
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthMax = 0.05;
Mesh.RemeshAlgorithm = 1; //(0) nosplit (1) automatic (2) split only with metis
Mesh 2;
Save StrCat(StrPrefix(General.FileName), ".vtk");
