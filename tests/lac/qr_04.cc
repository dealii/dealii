// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
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
 * Same qr_02.cc but test ImplicitQR::remove_column()
 */

/*
 * MWE in Python:

import numpy as np
from scipy import linalg
from math import sqrt

A = np.array([[0, 1, 2], [3, 4, 5], [6, 3, 4.5], [7,9, 8], [10,8,0]])
Q, R = linalg.qr(A)
Q = -Q
R = -R
# first element on diagonal should be
# sqrt(R[0,1]**2+R[1,1]**2) = 13.076696830622021
Q1, R1 = linalg.qr_delete(Q, R, 0, 1, 'col', False)
# thin QR
ncol = A.shape[1]-1
R1 = R1[:ncol,:]
Q1 = Q1[:, :ncol]
np.set_printoptions(precision=6)

>>> print A
[[  0.    1.    2. ]
 [  3.    4.    5. ]
 [  6.    3.    4.5]
 [  7.    9.    8. ]
 [ 10.    8.    0. ]]

>>> print R1
[[ 13.076697   8.22073 ]
 [  0.         6.757928]]

>>> print Q1
[[ 0.076472  0.202924]
 [ 0.305888  0.367773]
 [ 0.229416  0.38681 ]
 [ 0.688247  0.346572]
 [ 0.611775 -0.744198]]

# first Givens values for plane 0,1 are 0.949833 0.312758 13.0767
# so we expect the first column to be
# 0.949833*Q[:,0]+0.312758*Q[:,1]
*/


template <typename VectorType>
void
print(const ImplicitQR<VectorType> &qr, const unsigned int col_size)
{
  const LAPACKFullMatrix<double> &R = qr.get_R();
  std::vector<VectorType>         A(R.n());
  for (unsigned int j = 0; j < R.n(); ++j)
    {
      Vector<double> x(R.n());
      x    = 0;
      x[j] = 1.;
      A[j].reinit(col_size);
      qr.multiply_with_A(A[j], x);
    }

  const unsigned int size = qr.size();
  deallog << size << std::endl;
  deallog << "R:" << std::endl;
  R.print_formatted(deallog.get_file_stream(), 6, false, 10);
  deallog << "A:" << std::endl;

  for (unsigned int i = 0; i < col_size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      {
        deallog.get_file_stream() << std::setw(9) << A[j](i);
        if (j < size - 1)
          deallog.get_file_stream() << ' ';
        else
          deallog.get_file_stream() << std::endl;
      }
}

template <typename number>
void
test()
{
  using VectorType = Vector<number>;
  ImplicitQR<VectorType> qr;

  auto print_givens = [](const unsigned int           i,
                         const unsigned int           j,
                         const std::array<number, 3> &csr) {
    deallog.get_file_stream() << "Givens " << i << ' ' << j << ": " << csr[0]
                              << ' ' << csr[1] << ' ' << csr[2] << std::endl;
  };
  qr.connect_givens_slot(print_givens);


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

  // up to here ^^^^ exactly the same as in qr_03 test
  print(qr, v_size);

  // remove first column
  qr.remove_column();

  print(qr, v_size);
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  test<double>();
}
