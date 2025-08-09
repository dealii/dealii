/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------*/

r = 3; 

d1 = 0.1;
b1 = 0.3;
a2 = 0.5;
b2 = 0.7;
d2 = 1.0; 

Point(0) = { 0, 0, 0};
Point(1) = { d1, 0 , 0};
Point(2) = { d1, d1, 0};
Point(3) = { 0, d1, 0};

Point(4) = {b1, 0, 0};
Point(5) = {b1/Sqrt(2), b1/Sqrt(2), 0};
Point(6) = { 0, b1 , 0};

Point(7) = {a2, 0, 0};
Point(8) = {a2/Sqrt(2), a2/Sqrt(2), 0};
Point(9) = { 0, a2, 0};

Point(10) = {b2, 0, 0};
Point(11) = {b2/Sqrt(2), b2/Sqrt(2), 0};
Point(12) = { 0, b2, 0};

Point(13) = {d2, 0, 0};
Point(14) = {d2/Sqrt(2), d2/Sqrt(2), 0};
Point(15) = { 0, d2, 0};

Line(1) = {1, 4};
Line(2) = {4, 7};
Line(3) = {7, 10};
Line(4) = {10, 13};

Line(5) = {2, 5};
Line(6) = {5, 8};
Line(7) = {8, 11};
Line(8) = {11, 14};

Line(9) = {3, 6};
Line(10) = {6, 9};
Line(11) = {9, 12};
Line(12) = {12, 15};

Line(13) = {1, 2};
Line(14) = {2, 3};

Circle(15) = {4, 0, 5};
Circle(16) = {5, 0, 6};

Circle(17) = {7, 0, 8};
Circle(18) = {8, 0, 9};

Circle(19) = {10, 0, 11};
Circle(20) = {11, 0, 12};

Circle(21) = {13, 0, 14};
Circle(22) = {14, 0, 15};

Line Loop(1) = {1, 15, -5, -13};
Plane Surface(1) = {1};

Line Loop(2) = {2, 17, -6, -15};
Plane Surface(2) = {2};

Line Loop(3) = {3, 19, -7, -17};
Plane Surface(3) = {3};

Line Loop(4) = {4, 21, -8, -19};
Plane Surface(4) = {4};

Line Loop(5) = {5, 16, -9, -14};
Plane Surface(5) = {5};

Line Loop(6) = {6, 18, -10, -16};
Plane Surface(6) = {6};

Line Loop(7) = {7, 20, -11, -18};
Plane Surface(7) = {7};

Line Loop(8) = {8, 22, -12, -20};
Plane Surface(8) = {8};

Q1[] = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{1, 2, 3, 4, 5, 6, 7, 8};} };
Q2[] = Rotate { {0,0,1}, {0,0,0}, Pi}
{ Duplicata{ Surface{1,2,3,4,5,6,7,8,Q1[0],Q1[1],Q1[2],Q1[3],Q1[4],Q1[5],Q1[6],Q1[7]};} };

l = newl;  Line(l) = {0, 1};
l = newl;  Line(l) = {0, 3};
l = newl;  Line(l) = {0, 97};
l = newl;  Line(l) = {0, 228};

ll = newll; Line Loop(ll) = {13, 14, -141, 140};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {141, -27, -47, -142};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-143, 142, -66, -86};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-126, -140, 143, -106};
s = news; Plane Surface(s) = {ll};

Physical Surface(1) = {Surface{:}} ;

// The following sets the boundary ID. The parameter 
// Settings::outer_boundary_id in the program must equal the number in the
// brackets below, i.e., Settings::outer_boundary_id = 1. 
Physical Line(1) = {21, 22, 40, 60, 79, 99, 119, 139};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;
