// Gmsh project created on Fri Jul 26 13:43:43 2024
SetFactory("OpenCASCADE");
ms = 0.2;
Point(1) = {0, 3, 3, ms};
Point(2) = {0, 3, -3, ms};
Point(3) = {0, -3, -3, ms};
Point(4) = {0, -3, 3, ms};
Point(5) = {20, 3, 3, ms};
Point(6) = {20, 3, -3, ms};
Point(7) = {20, -3, -3, ms};
Point(8) = {20, -3, 3, ms};
Point(9) = {20, 3, 0, ms};
Point(10) = {20, -3, 0, ms};

//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Line(5) = {4, 8};
//+
Line(6) = {3, 7};
//+
Line(7) = {2, 6};
//+
Line(8) = {1, 5};
//+
Line(9) = {8, 5};
//+
Line(10) = {7, 6};
//+
Line(11) = {8, 10};
//+
Line(12) = {8, 10};
//+
Line(13) = {10, 9};
//+
Line(14) = {9, 5};
//+
Line(15) = {10, 7};
//+
Line(16) = {9, 6};
//+
Curve Loop(1) = {4, 5, 9, -8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 6, -15, -11, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 8, -14, 16, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, 7, -10, -6};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, 13, 14, -9};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {13, 16, -10, -15};
//+
Plane Surface(6) = {6};
//+


Field[1] = Box;
//+
Field[1].Thickness = 1.5;
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
Field[1].ZMax = -1.0;
//+
Field[1].ZMin = 1.0;
//+
Background Field = 0;





Physical Surface("fb", 17) = {2, 3};
//+
Curve Loop(7) = {4, 1, 2, 3};
//+
Plane Surface(7) = {7};
//+
Physical Surface("inlet", 18) = {7};
//+
Physical Surface("outlet", 19) = {5};
//+
Physical Surface("top", 20) = {1};
//+
Physical Surface("bottom", 21) = {4};
//+
Physical Surface("wall", 22) = {6};
//+
Surface Loop(1) = {7, 1, 2, 4, 3, 5, 6};
//+
Volume(1) = {1};
//+
Physical Volume("fluid", 23) = {1};


Mesh.SubdivisionAlgorithm = 0;
Mesh.MshFileVersion = 2.13;//+
