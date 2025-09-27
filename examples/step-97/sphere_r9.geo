/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------*/

r = 9;

d1 = 0.1;
a1 = 0.3;
b1 = 0.6;
a2 = 0.9;
b2 = 1.2;
d2 = 2.0;
d3 = 3.0;

a13 = a1/Sqrt(3);
b13 = b1/Sqrt(3);

a23 = a2/Sqrt(3);
b23 = b2/Sqrt(3);

d23 = d2/Sqrt(3);
d33 = d3/Sqrt(3);

Point(0) = { 0, 0, 0};

Point(1) = { d1, d1, d1};
Point(2) = { d1,-d1, d1};

Point(3) = { a13, a13, a13};
Point(4) = { a13,-a13, a13};

Point(5) = { b13, b13, b13};
Point(6) = { b13,-b13, b13};

Point(7) = { a23, a23, a23};
Point(8) = { a23,-a23, a23};

Point(9) = { b23, b23, b23};
Point(10) = { b23,-b23, b23};

Point(11) = { d23, d23, d23};
Point(12) = { d23,-d23, d23};

Point(13) = { d33, d33, d33};
Point(14) = { d33,-d33, d33};

Line(1) = {1, 3};
Line(2) = {3, 5};
Line(3) = {5, 7};
Line(4) = {7, 9};
Line(5) = {9, 11};

Line(7) = {2, 4};
Line(8) = {4, 6};
Line(9) = {6, 8};
Line(10) = {8, 10};
Line(11) = {10, 12};

Line(13) = {1, 2};

Circle(14) = {3, 0, 4};
Circle(15) = {5, 0, 6};
Circle(16) = {7, 0, 8};
Circle(17) = {9, 0, 10};

Circle(18) = {11, 0, 12};

Line(19) = {11,13};
Line(20) = {12,14};

Circle(21) = {13, 0, 14};

Line Loop(1) = {7, -14, -1, 13};
Plane Surface(1) = {1};

Line Loop(2) = {8, -15, -2, 14};
Plane Surface(2) = {2};

Line Loop(3) = {9, -16, -3, 15};
Plane Surface(3) = {3};

Line Loop(4) = {10, -17, -4, 16};
Plane Surface(4) = {4};

Line Loop(5) = {11, -18, -5, 17};
Plane Surface(5) = {5};

Line Loop(6) = {20, -21, -19, 18};
Plane Surface(6) = {6};

Q1[] = Symmetry {0, 1, 1, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {0, -1, 1, 0} {Duplicata {Surface{:};}};

ll = newll; Line Loop(ll) = {13, -26, 85, -55};
Plane Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {14, 24, -83, 53};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {15, 29, -88, 58};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {16, 34, -93, 63};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {17, 39, -98, 68};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {18, 44, -103, 73};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {21, 49, -108, 78};
Surface(ll) = {ll};

Surface Loop(1) = {1, 22, 81, 51, 109, 110};
Volume(1) = {1};

Surface Loop(2) = {2, 27, 86, 56, 110, 111};
Volume(2) = {2};

Surface Loop(3) = {3, 32, 91, 61, 111, 112};
Volume(3) = {3};

Surface Loop(4) = {4, 37, 96, 66, 112, 113};
Volume(4) = {4};

Surface Loop(5) = {5, 42, 101, 71, 113, 114};
Volume(5) = {5};

Surface Loop(6) = {6, 47, 106, 76, 114, 115};
Volume(6) = {6};

V1[] = Symmetry {1, 1, 0, 0} {Duplicata {Volume{:};}};
V2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Volume{:};}};
V3[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2, 3, 4, 5, 6};}};
V4[] = Symmetry {1, 0, -1, 0} {Duplicata {Volume{1, 2, 3, 4, 5, 6};}};

Surface Loop(7) = {109, 319, 505, 137, 687, 869};
Volume(7) = {7};

Physical Volume(1) = {Volume{:}};

Physical Surface(2) = {115, 479, 665, 297, 847, 1029};

Recombine Surface "*";
Recombine Volume "*";

Transfinite Volume "*";
Transfinite Surface "*";
Transfinite Line "*" = r;
