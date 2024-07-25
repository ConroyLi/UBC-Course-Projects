// Gmsh project created on Thu May 23 16:59:47 2024
SetFactory("OpenCASCADE");
Merge "domain_increase_xy.STEP";
//Merge "domain_tank.STEP";
//Merge "domain_increase_y.STEP";
//Merge "domain_reduce_y.STEP";
//Merge "ShipDomain_shift_x.STEP";
//Merge "ShipDomain_new.STEP";//+

S = 18;
Field[1] = Box;
//+
Field[1].Thickness = 1.5;
//+
Field[1].VIn = 0.1;
//+
Field[1].VOut = 0.85;
//+
Field[1].XMax = 21; //+S
//+
Field[1].XMin = -18;
//+
Field[1].YMax = 7; //5.72
//+
Field[1].YMin = -7;
//+
Field[1].ZMax = -0.8;
//+
Field[1].ZMin = 1.3;
//+
Background Field = 1;

//+
Physical Surface("inlet", 27) = {5};
//+
Physical Surface("top", 28) = {4};
//+
Physical Surface("bottom", 29) = {2};
//+
Physical Surface("outlet", 30) = {6};
//+
Physical Surface("ship", 31) = {10, 12, 8, 7, 11, 13, 9};
//+
Physical Surface("fb", 33) = {3, 1};
//+
Physical Volume("fluid", 32) = {1};

Mesh.SubdivisionAlgorithm = 0;
Mesh.MshFileVersion = 2.13;//+
