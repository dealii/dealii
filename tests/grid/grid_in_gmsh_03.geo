hmin=0.5;
Point (2) = {0, 0, 0 ,hmin };
Point (3) = {2, 0, 0 ,hmin };
Point (4) = {0, 2, 0 ,hmin };
Point (5) = {2, 2, 0 ,hmin };
Point (6) = {1, 0, 0 ,hmin };
Point (7) = {0, 1, 0 ,hmin };

Line (2) = {2, 7 };
Line (3) = {7, 4 };
Line (4) = {4, 5 };
Line (5) = {5, 3 };
Line (6) = {3, 6 };
Line (7) = {6, 2 };

Circle(8) = {6, 2, 7};
Curve Loop(1) = {7, 2, -8};
Plane Surface(1) = {1};

Curve Loop(2) = {8, 3, 4, 5, 6};
Plane Surface(2) = {2};

Physical Curve("1", 10) = {2, 7, 6, 5, 4, 3};
Physical Surface("1", 9) = {2, 1};
