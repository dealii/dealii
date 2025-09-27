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
 * Test QR::multiply_with_XY()
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
Q = -Q[:, :ncol]

d1 = np.array([1,3,5,7,9])
d2 = np.array([10,9,9])

np.set_printoptions(precision=6)

>>> A.dot(d2)
array([  27. ,  111. ,  127.5,  223. ,  172. ])

>>> A.T.dot(d1)
array([ 178. ,  163. ,   95.5])

>>> Q.dot(d2)
array([  3.342054,   7.775949,   5.279812,  12.737741,  -0.488702])

>>> Q.T.dot(d1)
array([ 12.779655,   1.043571,   0.071793])

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

  Vector<number> d1(5), d2(3), res1(3), res2(5);
  d1[0] = 1;
  d1[1] = 3;
  d1[2] = 5;
  d1[3] = 7;
  d1[4] = 9;

  d2[0] = 10;
  d2[1] = 9;
  d2[2] = 9;


  const unsigned int size = qr.size();
  deallog << size << std::endl;

  deallog << "A:" << std::endl;
  qr.multiply_with_A(res2, d2);
  res2.print(deallog.get_file_stream(), 6, false);

  deallog << "AT:" << std::endl;
  qr.multiply_with_AT(res1, d1);
  res1.print(deallog.get_file_stream(), 6, false);

  deallog << "Q:" << std::endl;
  qr.multiply_with_Q(res2, d2);
  res2.print(deallog.get_file_stream(), 6, false);

  deallog << "QT:" << std::endl;
  qr.multiply_with_QT(res1, d1);
  res1.print(deallog.get_file_stream(), 6, false);
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  test<double>();
}
