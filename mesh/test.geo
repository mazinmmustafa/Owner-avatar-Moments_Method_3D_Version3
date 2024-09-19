//+
Point(1) = {-1, -1, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, -1, 0, 1.0};
//+
Point(5) = {-1, 1, 0, 1.0};
//+

//+
Point(6) = {-0.5, -0.5, 0, 1.0};
//+
Point(7) = {-0.5, 0.5, 0, 1.0};
//+
Point(8) = {0.5, 0.5, 0, 1.0};
//+
Point(9) = {0.5, -0.5, 0, 1.0};
//+
Line(1) = {3, 5};
//+
Line(2) = {5, 1};
//+
Line(3) = {1, 4};
//+
Line(4) = {4, 3};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 9};
//+
Line(7) = {9, 6};
//+
Line(8) = {6, 7};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};