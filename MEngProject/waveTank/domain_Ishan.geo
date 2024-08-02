// Gmsh project created on Fri Apr 12 12:23:53 2024
fine_mesh = 1.0;
coarse_mesh = 1.0;

domain_length = 14;
domain_x_start = 0;
coarse_length = 5;

domain_width = 6;
domain_height = 6;

Point(1) = {domain_x_start, -domain_width/2, -domain_height/2, fine_mesh};
Point(2) = {domain_x_start, domain_width/2, -domain_height/2, fine_mesh};
//Point(3) = {domain_x_start + domain_length - coarse_length, domain_width/2, -domain_height/2, fine_mesh};
Point(4) = {domain_x_start + domain_length, domain_width/2, -domain_height/2, coarse_mesh};
Point(5) = {domain_x_start + domain_length, -domain_width/2, -domain_height/2, coarse_mesh};
//Point(6) = {domain_x_start + domain_length - coarse_length, -domain_width/2, -domain_height/2, fine_mesh};

Line(1) = {1,2};
//Line(2) = {2,3};
//Line(3) = {3,4};
Line(4) = {4,5};
//Line(5) = {5,6};
//Line(6) = {6,1};

Curve Loop(1) = {6, 1, 2, 3, 4, 5};
Plane Surface(1) = {1};

Extrude {0, 0, domain_height} {
  Surface{1};
}

// Block

block_length = 3;
block_depth = 0;
block_tr_length = 0.9;
block_start_pos = 1.78;

Point(25) = {block_start_pos, -block_tr_length/2, -block_tr_length/2 - block_depth, fine_mesh};
Point(26) = {block_start_pos, block_tr_length/2, -block_tr_length/2 - block_depth, fine_mesh};
Point(27) = {block_start_pos + block_length, block_tr_length/2, -block_tr_length/2 - block_depth, fine_mesh};
Point(28) = {block_start_pos + block_length, -block_tr_length/2, -block_tr_length/2 - block_depth, fine_mesh};

Line(33) = {25,26};
Line(34) = {26,27};
Line(35) = {27,28};
Line(36) = {28,25};

Curve Loop(2) = {33, 34, 35, 36};
Plane Surface(2) = {2};

Extrude {0, 0, block_tr_length} {
  Surface{2};
}

Delete{ Volume{1}; Volume{2}; }
Surface Loop(1) = {38, 17, 1, 21, 25, 29, 33, 37};
Surface Loop(2) = {47, 2, 51, 55, 59, 60};
Volume(1) = {1, 2};

Field[1] = Box;
Field[1].VIn = 0.1;
Field[1].VOut = 1;
Field[1].XMax = domain_x_start + domain_length - coarse_length;
Field[1].XMin = domain_x_start;
Field[1].YMax = domain_width/2;
Field[1].YMin = -domain_width/2;
Field[1].ZMax = 1;
Field[1].ZMin = -1;
Field[2] = Min;
Field[2].FieldsList = {1};
Background Field = 2;

Physical Surface("inlet") = {21};
Physical Surface("outlet") = {33};
Physical Surface("fb") = {25, 17, 29, 37};
Physical Surface("top") = {38};
Physical Surface("bottom") = {1};
Physical Surface("block") = {47, 51, 60, 55, 59, 2};
Physical Surface("block_face") = {47};
Physical Volume("fluid") = {1};

Mesh.SubdivisionAlgorithm = 0;
Mesh.MshFileVersion = 2.13;
