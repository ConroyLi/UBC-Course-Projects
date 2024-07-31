// Gmsh project created on Thu May 23 16:59:47 2024
SetFactory("OpenCASCADE");

//+
Point(1) = {0, 3, 3, 1.0};
//+
Point(2) = {0, 3, -3, 1.0};
//+
Point(3) = {0, -3, -3, 1.0};
//+
Point(4) = {0, -3, 3, 1.0};
//+
Point(5) = {20, -3, 3, 1.0};
//+
Point(6) = {20, -3, -3, 1.0};
//+
Point(7) = {20, 3, -3, 1.0};
//+
Point(8) = {20, 3, 3, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {1, 8};
//+
Line(7) = {2, 7};
//+
Line(8) = {3, 6};
//+
Line(9) = {5, 6};
//+
Line(10) = {6, 7};
//+
Line(11) = {7, 8};
//+
Line(12) = {8, 5};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 7, 11, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 6, 12, -5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 10, -7, 3};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8, -9, -5, -4};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {11, 12, 9, 10};
//+
Plane Surface(6) = {6};

Field[1] = Box;
//+
Field[1].Thickness = 8;
//+
Field[1].VIn = 0.1;
//+
Field[1].VOut = 1;
//+
Field[1].XMax = 15; 
//+
Field[1].XMin = 0;
//+
Field[1].YMax = 3; 
//+
Field[1].YMin = -3;
//+
Field[1].ZMax = -0.8;
//+
Field[1].ZMin = 1.3;
//+
Background Field = 1;

//+
Point(9) = {1.78, 0.7,  0.7, 1.0};
//+
Point(10) = {1.78,  0.7, - 0.7, 1.0};
//+
Point(11) = {1.78, - 0.7, - 0.7, 1.0};
//+
Point(12) = {1.78, - 0.7,  0.7, 1.0};
//+
Point(13) = {4.78, - 0.7,  0.7, 1.0};
//+
Point(14) = {4.78, - 0.7, - 0.7, 1.0};
//+
Point(15) = {4.78,  0.7, - 0.7, 1.0};
//+
Point(16) = {4.78,  0.7,  0.7, 1.0};
//+
Line(13) = {10, 9};
//+
Line(14) = {9, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {11, 10};
//+
Line(17) = {9, 16};
//+
Line(18) = {12, 13};
//+
Line(19) = {10, 15};
//+
Line(20) = {11, 14};
//+
Line(21) = {16, 15};
//+
Line(22) = {15, 14};
//+
Line(23) = {14, 13};
//+
Line(24) = {13, 16};
//+
Curve Loop(7) = {15, 16, 13, 14};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {14, 18, 24, -17};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {15, 20, 23, -18};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {17, 21, -19, 13};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {19, 22, -20, 16};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {23, 24, 21, 22};
//+
Plane Surface(12) = {12};
//+
Physical Surface("inlet", 25) = {1};
//+
Physical Surface("top", 26) = {3};
//+
Physical Surface("bottom", 27) = {4};
//+
Physical Surface("outlet", 28) = {6};
//+
Physical Surface("fb", 29) = {5, 2};
//+
Physical Surface("ship", 30) = {7, 8, 10, 9, 11, 12};
//+
Surface Loop(1) = {5, 4, 6, 2, 1, 3};
//+
Surface Loop(2) = {9, 7, 11, 10, 8, 12};
//+
Volume(1) = {1, 2};
//+
Physical Volume("fluid", 31) = {1};
