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


// Same as 'tensor_product_matrix_06.cc' except that we replaced the scalar data
// type 'float' by the vectorized data type 'VectorizedArray<float>'.
// Note, all lanes compute the same.

#include <deal.II/base/vectorization.h>

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
  Table<2, VectorizedArray<float>> init_mass(size, size);
  Table<2, VectorizedArray<float>> init_laplace(size, size);
  for (unsigned int i = 0; i < size; ++i)
    {
      init_mass(i, i) = 2. / 3.;
      if (i > 0)
        init_mass(i, i - 1) = 1. / 6.;
      if (i < size - 1)
        init_mass(i, i + 1) = 1. / 6.;
      init_laplace(i, i) = 2.;
      if (i > 0)
        init_laplace(i, i - 1) = -1.;
      if (i < size - 1)
        init_laplace(i, i + 1) = -1.;
    }
  std::array<Table<2, VectorizedArray<float>>, dim> mass;
  std::array<Table<2, VectorizedArray<float>>, dim> laplace;
  for (unsigned int dir = 0; dir < dim; ++dir)
    {
      for (unsigned int i = 0; i < size; ++i)
        {
          init_mass(i, i) *= make_vectorized_array<float>(4. / 3.);
          init_laplace(i, i) *= make_vectorized_array<float>(5. / 4.);
        }
      mass[dir]    = init_mass;
      laplace[dir] = init_laplace;
    }
  TensorProductMatrixSymmetricSum<dim, VectorizedArray<float>> mat;
  mat.reinit(mass, laplace);

  Vector<float> w1(mat.m()), w2(mat.m());
  for (unsigned int i = 0; i < w1.size(); ++i)
    w1[i] = (2 * i + 1) % 23;

  auto convert_to_vectorized = [](const Vector<float> &                  in,
                                  AlignedVector<VectorizedArray<float>> &out) {
    std::transform(in.begin(), in.end(), out.begin(), [](const float &val) {
      return make_vectorized_array(val);
    });
  };
  AlignedVector<VectorizedArray<float>> v1(w1.size()), v2(w1.size()),
    v3(w1.size());
  convert_to_vectorized(w1, v1);

  constexpr unsigned int macro_size = VectorizedArray<float>::size();
  Vector<float>          vec_flat(v1.size() * macro_size);
  std::array<unsigned int, macro_size> offsets;
  for (unsigned int i = 0; i < macro_size; ++i)
    offsets[i] = v1.size() * i;
  auto subtract_and_assign =
    [](AlignedVector<VectorizedArray<float>> &      lhs,
       const AlignedVector<VectorizedArray<float>> &rhs) {
      std::transform(lhs.begin(),
                     lhs.end(),
                     rhs.begin(),
                     lhs.begin(),
                     [](const VectorizedArray<float> lval,
                        const VectorizedArray<float> rval) {
                       return lval - rval;
                     });
    };

  const ArrayView<VectorizedArray<float>> view1(v1.begin(), v1.size());
  const ArrayView<VectorizedArray<float>> view2(v2.begin(), v2.size());
  const ArrayView<VectorizedArray<float>> view3(v3.begin(), v3.size());
  mat.vmult(view2, view1);
  mat.apply_inverse(view3, view2);
  subtract_and_assign(v3, v1);
  vectorized_transpose_and_store(
    false, mat.m(), v3.begin(), offsets.data(), vec_flat.begin());
  double norm = vec_flat.linfty_norm();
  deallog << "Verification of vmult and inverse: " << (norm < 1e-3 ? 0. : norm)
          << std::endl;

  FullMatrix<float> full(v1.size(), v1.size());
  full = 0.;
  for (unsigned int dir = 0; dir < dim; ++dir)
    for (unsigned int i = 0, c = 0; i < (dim > 2 ? size : 1); ++i)
      for (unsigned int j = 0; j < (dim > 1 ? size : 1); ++j)
        for (unsigned int k = 0; k < size; ++k, ++c)
          for (unsigned int ii = 0, cc = 0; ii < (dim > 2 ? size : 1); ++ii)
            for (unsigned int jj = 0; jj < (dim > 1 ? size : 1); ++jj)
              for (unsigned int kk = 0; kk < size; ++kk, ++cc)
                if (dim == 1)
                  full(c, cc) = laplace[0](k, kk)[0];
                else if (dim == 2)
                  full(c, cc) = laplace[1](j, jj)[0] * mass[0](k, kk)[0] +
                                mass[1](j, jj)[0] * laplace[0](k, kk)[0];
                else if (dim == 3)
                  full(c, cc) = laplace[2](i, ii)[0] * mass[1](j, jj)[0] *
                                  mass[0](k, kk)[0] +
                                mass[2](i, ii)[0] *
                                  (laplace[1](j, jj)[0] * mass[0](k, kk)[0] +
                                   mass[1](j, jj)[0] * laplace[0](k, kk)[0]);
  full.vmult(w2, w1);
  convert_to_vectorized(w2, v3);
  subtract_and_assign(v3, v2);
  vectorized_transpose_and_store(
    false, mat.m(), v3.begin(), offsets.data(), vec_flat.begin());
  norm = vec_flat.linfty_norm();
  deallog << "Verification of vmult: " << (norm < 1e-4 ? 0. : norm)
          << std::endl;

  full.gauss_jordan();
  full.vmult(w2, w1);
  mat.apply_inverse(view2, view1);
  convert_to_vectorized(w2, v3);
  subtract_and_assign(v3, v2);
  vectorized_transpose_and_store(
    false, mat.m(), v3.begin(), offsets.data(), vec_flat.begin());
  norm = vec_flat.linfty_norm();
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
