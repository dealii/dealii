// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/qr.h>

#include "../tests.h"

/*
 * Test QR::append_column()
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
  std::vector<VectorType>         Q(qr.size());
  for (unsigned int j = 0; j < qr.size(); ++j)
    {
      Vector<double> x(qr.size());
      x    = 0;
      x[j] = 1.;
      Q[j].reinit(v_size);
      qr.multiply_with_Q(Q[j], x);
    }

  const unsigned int size = qr.size();
  deallog << size << std::endl;
  deallog << "R:" << std::endl;
  R.print_formatted(deallog.get_file_stream(), 6, false, 10);
  deallog << "Q:" << std::endl;

  for (unsigned int i = 0; i < v_size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      {
        deallog.get_file_stream() << std::setw(9) << Q[j](i);
        if (j < size - 1)
          deallog.get_file_stream() << " ";
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
