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

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/qr.h>

#include "../tests.h"

/*
 * Same as qr_08, but test ImplicitQR class, namely
 * use signal for givens rotation to update factorizations Q^T H Q and Q^T u
 * when the first column is removed from the subspace.
 * We use QR class to build Q and project matrix and vector on the subspace,
 * but attach signal to ImplicitQR.
 */

/*
 * MWE in Python:

import numpy as np
from scipy import linalg
from math import sqrt

# vector subspace of R^6:
A = np.array([[0, 1, 2, 7], [3, 4, 5, 1], [6, 3, 4.5, 2.2], [7,9, 8, 0],
[10,8,0, 1.1], [7, 1, 9, 5]]) Q, R = linalg.qr(A) Q = -Q R = -R ncol =
A.shape[1] # thin QR: R1 = R[:ncol,:] Q1 = Q[:, :ncol]

# 1D Laplace matrix and some vector to be represented in the subspace
M =
np.array([[2,-1,0,0,0,0],[-1,2,-1,0,0,0],[0,-1,2,-1,0,0],[0,0,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
V = np.array([1,2,2,1,5,9])

# projection
H1 = Q1.T.dot(M.dot(Q1))
V1 = Q1.T.dot(V)

# now delete the first column
Q2, R2 = linalg.qr_delete(Q, R, 0, 1, 'col', False)
R2 = R2[:ncol-1,:]
Q2 = Q2[:, :ncol-1]

# updated projections:
H2 = Q2.T.dot(M.dot(Q2))
V2 = Q2.T.dot(V)

np.set_printoptions(precision=6)

>>> print A
[[  0.    1.    2.    7. ]
 [  3.    4.    5.    1. ]
 [  6.    3.    4.5   2.2]
 [  7.    9.    8.    0. ]
 [ 10.    8.    0.    1.1]
 [  7.    1.    9.    5. ]]

>>> print Q1
[[-0.        0.160817  0.221587 -0.955862]
 [ 0.19245   0.285897  0.335821  0.033159]
 [ 0.3849   -0.232291  0.045103  0.005283]
 [ 0.44905   0.613487  0.388792  0.239359]
 [ 0.6415    0.095299 -0.70425  -0.16568 ]
 [ 0.44905  -0.673048  0.434697 -0.021413]]

>>> print R1
[[ 15.588457  11.547005  10.328155   3.990132]
 [ -0.         6.218253  -0.443735  -2.359839]
 [ -0.        -0.         9.347851   3.384966]
 [ -0.        -0.        -0.        -6.935562]]

>>> print H1
[[ 0.353909 -0.275486 -0.246262  0.084659]
 [-0.275486  2.337236 -0.110024  0.295859]
 [-0.246262 -0.110024  2.945693  0.587455]
 [ 0.084659  0.295859  0.587455  2.132729]]

>>> print V1
[ 8.852704 -4.699427  1.76325  -1.660732]

>>> print Q2
[[ 0.076249  0.123157  0.867562]
 [ 0.304997  0.213292 -0.108274]
 [ 0.228748  0.229803  0.073457]
 [ 0.686244  0.177292 -0.3507  ]
 [ 0.609994 -0.504539  0.294887]
 [ 0.076249  0.774943  0.142366]]

>>> print H2
[[ 0.569767 -0.567113 -0.321939]
 [-0.567113  2.728827 -0.480842]
 [-0.321939 -0.480842  2.378169]]

>>> print V2
[ 5.566198  5.638437  3.202953]
*/


template <typename VectorType>
void
print(const ImplicitQR<VectorType> &qr)
{
  const LAPACKFullMatrix<double> &R = qr.get_R();

  const unsigned int size = qr.size();
  deallog << size << std::endl;
  deallog << "R:" << std::endl;
  R.print_formatted(deallog.get_file_stream(), 6, false, 10);
}

template <typename number>
void
test()
{
  using VectorType = Vector<number>;
  QR<VectorType>         qr;
  ImplicitQR<VectorType> pqr;

  const unsigned int v_size = 6;
  VectorType         v(v_size);
  v[0] = 0;
  v[1] = 3;
  v[2] = 6;
  v[3] = 7;
  v[4] = 10.;
  v[5] = 7.;

  qr.append_column(v);
  pqr.append_column(v);

  v[0] = 1;
  v[1] = 4;
  v[2] = 3;
  v[3] = 9;
  v[4] = 8.;
  v[5] = 1.;

  qr.append_column(v);
  pqr.append_column(v);

  v[0] = 2;
  v[1] = 5;
  v[2] = 4.5;
  v[3] = 8;
  v[4] = 0;
  v[5] = 9.;

  qr.append_column(v);
  pqr.append_column(v);

  v[0] = 7;
  v[1] = 1;
  v[2] = 2.2;
  v[3] = 0;
  v[4] = 1.1;
  v[5] = 5.;

  qr.append_column(v);
  pqr.append_column(v);

  print(pqr);

  // matrix to be projected:
  FullMatrix<number> M(v_size);
  M       = 0.;
  M(0, 0) = 2.;
  M(0, 1) = -1;
  M(1, 0) = -1;
  M(1, 1) = 2;
  M(1, 2) = -1;
  M(2, 1) = -1;
  M(2, 2) = 2;
  M(2, 3) = -1;
  M(3, 2) = -1;
  M(3, 3) = 2;
  M(3, 4) = -1;
  M(4, 3) = -1;
  M(4, 4) = 2;
  M(4, 5) = -1;
  M(5, 4) = -1;
  M(5, 5) = 2;

  // vector to be projected:
  Vector<number> V(v_size);
  V[0] = 1;
  V[1] = 2;
  V[2] = 2;
  V[3] = 1;
  V[4] = 5;
  V[5] = 9;

  // reduced matrix
  // H_{ij} = Q(k,i) M(k,l) Q(l,j)
  std::vector<VectorType> Q(qr.size());
  for (unsigned int j = 0; j < qr.size(); ++j)
    {
      Vector<double> x(qr.size());
      x    = 0;
      x[j] = 1.;
      Q[j].reinit(v_size);
      qr.multiply_with_Q(Q[j], x);
    }

  LAPACKFullMatrix<number> H1(4);
  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int j = 0; j < 4; ++j)
      {
        const VectorType &qi  = Q[i];
        const VectorType &qj  = Q[j];
        number            res = 0.;
        for (unsigned int k = 0; k < v_size; ++k)
          for (unsigned int l = 0; l < v_size; ++l)
            res += qi[k] * M(k, l) * qj[l];

        H1(i, j) = res;
      }

  deallog.get_file_stream() << std::endl;
  deallog << "H1:" << std::endl;
  H1.print_formatted(deallog.get_file_stream(), 6, false, 10);

  Vector<number> V1(4);
  qr.multiply_with_QT(V1, V);
  deallog << "V1:" << std::endl;
  V1.print(deallog.get_file_stream(), 6, false);

  // remove first column
  auto update_factorization = [&](const unsigned int           i,
                                  const unsigned int           k,
                                  const std::array<number, 3> &csr) {
    H1.apply_givens_rotation(csr, i, k, true);
    H1.apply_givens_rotation(csr, i, k, false);
    V1.apply_givens_rotation(csr, i, k);
  };

  pqr.connect_givens_slot(update_factorization);
  pqr.remove_column();

  // remove one row (and column)
  Vector<number>           V2(3);
  LAPACKFullMatrix<number> H2(3);
  for (unsigned int i = 0; i < 3; ++i)
    {
      V2(i) = V1(i);
      for (unsigned int j = 0; j < 3; ++j)
        H2(i, j) = H1(i, j);
    }


  deallog.get_file_stream() << std::endl;
  deallog << "H2:" << std::endl;
  H2.print_formatted(deallog.get_file_stream(), 6, false, 10);

  deallog << "V2:" << std::endl;
  V2.print(deallog.get_file_stream(), 6, false);
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  test<double>();
}
