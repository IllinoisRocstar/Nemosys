//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.2, 1.2, 0, 1.0};
//+
Recursive Delete {
  Point{5}; 
}
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {1, 3};
//+
Line(6) = {4, 2};
//+
Curve Loop(1) = {1, 2, -5};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, 3};
//+
Surface(2) = {2};
//+
Recursive Delete {
  Curve{6}; 
}
