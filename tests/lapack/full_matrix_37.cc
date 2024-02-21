// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Test application of Givens rotations to matrix and vector
 */

#include <deal.II/base/logstream.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

/*
 * MWE in Python:

import numpy as np
from scipy import linalg

A = np.array([[1.1, 1, 2, 7], [3, 4, 5, 1], [6, 3, 4.5, 2.2], [7,9, 8, 2.3]])
V = np.array([3,5,2,1])
c = 0.949833
s = 0.312758

i = 1
k = 3
G = np.array([[1,0,0,0], [0,c,0,s], [0,0,1,0], [0,-s,0,c]])

# from left
A1 = G.dot(A)

# from right
A2 = A.dot(G.T)

# vector
V1 = G.dot(V)

np.set_printoptions(precision=6)

>>> print A
[[ 1.1  1.   2.   7. ]
 [ 3.   4.   5.   1. ]
 [ 6.   3.   4.5  2.2]
 [ 7.   9.   8.   2.3]]

>>> print A1
[[ 1.1       1.        2.        7.      ]
 [ 5.038805  6.614154  7.251229  1.669176]
 [ 6.        3.        4.5       2.2     ]
 [ 5.710557  7.297465  6.034874  1.871858]]

 >>> print A2
[[ 1.1       3.139139  2.        6.336073]
 [ 3.        4.11209   5.       -0.301199]
 [ 6.        3.537567  4.5       1.151359]
 [ 7.        9.26784   8.       -0.630206]]

>>> print V1
[ 3.        5.061923  2.       -0.613957]

*/


template <typename NumberType>
void
test()
{
  const std::array<NumberType, 3> csr = {{0.949833, 0.312758, 13.0767}};
  LAPACKFullMatrix<NumberType>    A(4);
  Vector<NumberType>              V(4);
  const unsigned int              i = 1;
  const unsigned int              k = 3;

  A(0, 0) = 1.1;
  A(0, 1) = 1;
  A(0, 2) = 2;
  A(0, 3) = 7;
  A(1, 0) = 3;
  A(1, 1) = 4;
  A(1, 2) = 5;
  A(1, 3) = 1;
  A(2, 0) = 6;
  A(2, 1) = 3;
  A(2, 2) = 4.5;
  A(2, 3) = 2.2;
  A(3, 0) = 7;
  A(3, 1) = 9;
  A(3, 2) = 8;
  A(3, 3) = 2.3;

  V(0) = 3;
  V(1) = 5;
  V(2) = 2;
  V(3) = 1;

  LAPACKFullMatrix<NumberType> A1(A), A2(A);
  Vector<NumberType>           V1(V);

  A1.apply_givens_rotation(csr, i, k, true);
  A2.apply_givens_rotation(csr, i, k, false);
  V1.apply_givens_rotation(csr, i, k);

  deallog << "A1:" << std::endl;
  A1.print_formatted(deallog.get_file_stream(), 6, false, 10);

  deallog << "A2:" << std::endl;
  A2.print_formatted(deallog.get_file_stream(), 6, false, 10);

  deallog << "V1:" << std::endl;
  V1.print(deallog.get_file_stream(), 6, false);
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  test<double>();
}
