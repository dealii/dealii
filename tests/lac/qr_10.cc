// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
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
 * Test ImplicitQR::append_column() with linearly dependent columns.
 * The output should be the same as qr_03 test plus extra output due to
 * the rejected column.
 */

/*
 * MWE in Python for standard QR with 3 linearly independent columns. We will
add the fourth one:

import numpy as np
from scipy import linalg

A = np.array([[0, 1, 2], [3, 4, 5], [6, 3, 4.5], [7,9, 8], [10,8,0]])
Q, R = linalg.qr(A)
# thin QR
ncol = A.shape[1]
R = -R[:ncol,:]
Q = -Q[:, :ncol]
np.set_printoptions(precision=6)

print Q
print R

>>> print A
[[  0.    1.    2. ]
 [  3.    4.    5. ]
 [  6.    3.    4.5]
 [  7.    9.    8. ]
 [ 10.    8.    0. ]]

>>> print R
[[ 13.928388  12.420676   7.03599 ]
 [ -0.         4.089842   4.916632]
 [ -0.        -0.         6.290594]]
>>> print Q
[[-0.        0.244508  0.126831]
 [ 0.215387  0.32391   0.300765]
 [ 0.430775 -0.57472   0.682727]
 [ 0.502571  0.674288  0.182604]
 [ 0.717958 -0.224343 -0.627689]]

*/


template <typename number>
void
test()
{
  using VectorType = Vector<number>;
  ImplicitQR<VectorType> qr;

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

  // Now try adding the second column again
  v[0] = 1;
  v[1] = 4;
  v[2] = 3;
  v[3] = 9;
  v[4] = 8.;

  const auto fun = [&](const Vector<number> &u,
                       const number         &rho2,
                       const number         &col_l2_norm2) {
    deallog << "Linearly dependent column" << std::endl;
    deallog << "u:" << std::endl;
    u.print(deallog.get_file_stream());
    deallog << "rho2:" << std::endl;
    deallog.get_file_stream() << rho2 << std::endl;
    deallog << "col2:" << std::endl;
    deallog.get_file_stream() << col_l2_norm2 << std::endl;
    const bool reject =
      rho2 < dealii::Utilities::fixed_power<2>(1e-4) * col_l2_norm2;
    deallog << "reject:" << std::endl;
    deallog.get_file_stream() << reject << std::endl;
    return !reject;
  };

  qr.connect_append_column_slot(fun);
  const bool res = qr.append_column(v);
  deallog << "Success: " << res << std::endl;

  const LAPACKFullMatrix<double> &R = qr.get_R();
  std::vector<VectorType>         A(qr.size());
  for (unsigned int j = 0; j < qr.size(); ++j)
    {
      Vector<double> x(qr.size());
      x    = 0;
      x[j] = 1.;
      A[j].reinit(v_size);
      qr.multiply_with_A(A[j], x);
    }

  const unsigned int size = qr.size();
  deallog << size << std::endl;
  deallog << "R:" << std::endl;
  R.print_formatted(deallog.get_file_stream(), 6, false, 10);
  deallog << "A:" << std::endl;

  for (unsigned int i = 0; i < v_size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      {
        deallog.get_file_stream() << std::setw(9) << A[j](i);
        if (j < size - 1)
          deallog.get_file_stream() << ' ';
        else
          deallog.get_file_stream() << std::endl;
      }
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  test<double>();
}
