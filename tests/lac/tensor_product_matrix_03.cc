// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test templated path of TensorProductMatrix with float instead of double

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/tensor_product_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"


template <int dim, int size>
void
do_test()
{
  deallog << "Testing dim=" << dim << ", degree=" << size << std::endl;
  FullMatrix<float> mass(size, size);
  FullMatrix<float> laplace(size, size);
  for (unsigned int i = 0; i < size; ++i)
    {
      mass(i, i) = 2. / 3.;
      if (i > 0)
        mass(i, i - 1) = 1. / 6.;
      if (i < size - 1)
        mass(i, i + 1) = 1. / 6.;
      laplace(i, i) = 2.;
      if (i > 0)
        laplace(i, i - 1) = -1.;
      if (i < size - 1)
        laplace(i, i + 1) = -1.;
    }
  TensorProductMatrixSymmetricSum<dim, float, size> mat;
  mat.reinit(mass, laplace);
  Vector<float> v1(mat.m()), v2(mat.m()), v3(mat.m());
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = (2 * i + 1) % 23;

  const ArrayView<float> view1(v1.begin(), v1.size());
  const ArrayView<float> view2(v2.begin(), v2.size());
  const ArrayView<float> view3(v3.begin(), v3.size());
  mat.vmult(view2, view1);
  mat.apply_inverse(view3, view2);
  v3 -= v1;
  // add tolerance to account for different BLAS/LAPACK combinations. Float is
  // too inaccurate so numdiff does not work...
  double norm = v3.linfty_norm();
  deallog << "Verification of vmult and inverse: " << (norm < 1e-3 ? 0. : norm)
          << std::endl;

  FullMatrix<float> full(v1.size(), v1.size());
  for (unsigned int i = 0, c = 0; i < (dim > 2 ? size : 1); ++i)
    for (unsigned int j = 0; j < (dim > 1 ? size : 1); ++j)
      for (unsigned int k = 0; k < size; ++k, ++c)
        for (unsigned int ii = 0, cc = 0; ii < (dim > 2 ? size : 1); ++ii)
          for (unsigned int jj = 0; jj < (dim > 1 ? size : 1); ++jj)
            for (unsigned int kk = 0; kk < size; ++kk, ++cc)
              if (dim == 1)
                full(c, cc) = laplace(k, kk);
              else if (dim == 2)
                full(c, cc) =
                  laplace(k, kk) * mass(j, jj) + laplace(j, jj) * mass(k, kk);
              else if (dim == 3)
                full(c, cc) = (laplace(k, kk) * mass(j, jj) +
                               laplace(j, jj) * mass(k, kk)) *
                                mass(i, ii) +
                              laplace(i, ii) * mass(j, jj) * mass(k, kk);
  full.vmult(v3, v1);
  v3 -= v2;

  norm = v3.linfty_norm();
  deallog << "Verification of vmult: " << (norm < 1e-4 ? 0. : norm)
          << std::endl;

  full.gauss_jordan();
  full.vmult(v3, v1);
  mat.apply_inverse(view2, view1);
  v3 -= v2;
  norm = v3.linfty_norm();
  deallog << "Verification of inverse: " << (norm < 5e-3 ? 0. : norm)
          << std::endl;
}


int
main()
{
  initlog();

  do_test<1, 1>();
  do_test<1, 2>();
  do_test<1, 5>();
  do_test<2, 1>();
  do_test<2, 2>();
  do_test<2, 5>();
  do_test<2, 11>();
  do_test<3, 1>();
  do_test<3, 2>();
  do_test<3, 3>();
  do_test<3, 7>();

  return 0;
}
