// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/qr.h>

#include "../tests.h"

/*
 * Test BaseQR::solve()
 * QR is the same as in qr.cc
 */

/*
 * MWE in Python:

import numpy as np
from scipy import linalg

A = np.array([[0, 1, 2], [3, 4, 5], [6, 3, 4.5], [7,9, 8], [10,8,0]])
Q, R = linalg.qr(A)
# thin QR
ncol = A.shape[1]
R = -R[:ncol,:]
np.set_printoptions(precision=6)
d = np.array([1,3,5])
y = linalg.solve_triangular(R, d)
y2 = linalg.solve_triangular(R, d, trans='T')

>>> print R
[[ 13.928388  12.420676   7.03599 ]
 [ -0.         4.089842   4.916632]
 [ -0.        -0.         6.290594]]

>>> print y
[-0.131756 -0.221995  0.794838]

>>> print y2
[ 0.071796  0.515484  0.31164 ]
*/

template <typename number>
void
test()
{
  using VectorType = Vector<number>;
  QR<VectorType> qr;

  const unsigned int v_size = 5;
  VectorType         v(v_size);
  v[0] = 0;
  v[1] = 3;
  v[2] = 6;
  v[3] = 7;
  v[4] = 10.;

  qr.append_column(v);

  v[0] = 1;
  v[1] = 4;
  v[2] = 3;
  v[3] = 9;
  v[4] = 8.;

  qr.append_column(v);

  v[0] = 2;
  v[1] = 5;
  v[2] = 4.5;
  v[3] = 8;
  v[4] = 0;

  qr.append_column(v);


  const LAPACKFullMatrix<double> &R = qr.get_R();

  const unsigned int size = qr.size();
  deallog << size << std::endl;
  deallog << "R:" << std::endl;
  R.print_formatted(deallog.get_file_stream(), 6, false, 10);

  Vector<double> x(size), y(size);
  y[0] = 1.;
  y[1] = 3.;
  y[2] = 5.;

  qr.solve(x, y, false);
  x.print(deallog.get_file_stream(), 6, false);

  qr.solve(x, y, true);
  x.print(deallog.get_file_stream(), 6, false);
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  test<double>();
}
