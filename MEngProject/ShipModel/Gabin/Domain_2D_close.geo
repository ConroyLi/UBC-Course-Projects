SetFactory("OpenCASCADE");
Merge "Domain_2D_close.step";

Transfinite Curve {23, 22} = 10 Using Progression 1;
Transfinite Curve {55, 81} = 10 Using Progression 1;
Transfinite Curve {138, 130} = 10 Using Progression 1;
Transfinite Curve {165, 142} = 10 Using Progression 1;


Surface Loop(3) = {77, 76, 75, 78, 79, 74};
Surface Loop(4) = {4, 15, 26, 16, 5, 21, 22, 41, 34, 73, 64, 25, 42, 69, 61, 68, 71, 63, 62, 50, 72, 67, 66, 45, 40, 51, 65, 53, 48, 13, 37, 32, 28, 44, 52, 54, 38, 36, 47, 39, 46, 43, 27, 31, 7, 6, 14, 30, 23, 29, 33, 9, 11, 35, 8, 2, 1, 3, 10, 24, 12, 70, 49, 55, 56, 57, 58, 59, 60, 20, 19, 18, 17};
Volume(3) = {3, 4};
Delete {
  Volume{2}; Volume{1}; 
}


Physical Surface("Inlet", 194) = {74};
Physical Surface("Top", 195) = {77};
Physical Surface("Bottom", 196) = {75};
Physical Surface("Fb", 197) = {76, 78};
Physical Surface("Outlet", 198) = {79};
Physical Surface("Ship", 199) = {1, 3, 4, 73, 34, 41, 69, 37, 38, 67, 72, 10, 66, 36, 70, 65, 71, 33, 8, 2, 45, 40, 35, 51, 49, 29, 52, 48, 39, 62, 50, 64, 9, 11, 30, 24, 14, 46, 43, 31, 27, 23, 7, 6, 47, 53, 12, 26, 15, 5, 63, 68, 13, 28, 44, 16, 17, 18, 19, 20, 21, 22, 32, 55, 56, 57, 58, 59, 60, 61, 54, 42, 25};
Physical Volume("Fluid", 200) = {3};


Field[1] = Box;
Field[1].VIn = 1;
Field[1].VOut = 10;
Field[1].XMax = 75;
Field[1].XMin = -125;
Field[1].YMax = 50;
Field[1].YMin = -50;
Field[1].ZMax = 25;
Field[1].ZMin = -25;
Background Field = 1;


Mesh.MshFileVersion = 2.13;
