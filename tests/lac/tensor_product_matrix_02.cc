// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// Test templated path of TensorProductMatrix

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
  FullMatrix<double> mass(size, size);
  FullMatrix<double> laplace(size, size);
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
  TensorProductMatrixSymmetricSum<dim, double, size> mat;
  mat.reinit(mass, laplace);
  Vector<double> v1(mat.m()), v2(mat.m()), v3(mat.m());
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = (2 * i + 1) % 23;

  const ArrayView<double> view1(v1.begin(), v1.size());
  const ArrayView<double> view2(v2.begin(), v2.size());
  const ArrayView<double> view3(v3.begin(), v3.size());
  mat.vmult(view2, view1);
  mat.apply_inverse(view3, view2);
  v3 -= v1;
  deallog << "Verification of vmult and inverse: " << v3.linfty_norm()
          << std::endl;

  FullMatrix<double> full(v1.size(), v1.size());
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
  deallog << "Verification of vmult: " << v3.linfty_norm() << std::endl;

  full.gauss_jordan();
  full.vmult(v3, v1);
  mat.apply_inverse(view2, view1);
  v3 -= v2;
  deallog << "Verification of inverse: " << v3.linfty_norm() << std::endl;
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
