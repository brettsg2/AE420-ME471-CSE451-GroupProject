// Gmsh project created on Sun Dec 02 19:59:05 2018
Merge "Flange.igs";
Merge "gasket.igs";//+
Curve Loop(1) = {3, 2, -1, 8, -7, -6, -5, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 10, -11, 12};
//+
Plane Surface(2) = {2};
//+
Physical Curve("FlangeID") = {3};
//+
Physical Curve("Gasket_TopFace") = {9};
//+

