Point(1) = {0, 0, 0, 1.0};
Extrude {1, 1, 0} {
  Point{1}; Layers{10}; Recombine;
}
Extrude {-1, 1, 0} {
  Line{1}; Layers{10}; Recombine;
}
Extrude {0, 0, 1.4142135623730951e+00} {
  Surface{5}; Layers{10}; Recombine;
}
Physical Surface(1) = {27, 5, 26, 14, 18, 22};
Physical Volume(2) = {1};
