Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.25};
Point(3) = {1, 1, 0, 0.5};
Point(4) = {0, 1, 0, 0.25};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Recombine Surface(6);
Physical Line(1) = {1, 2, 3, 4};
Physical Surface(2) = {6};
