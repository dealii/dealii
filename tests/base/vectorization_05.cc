// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test transpose operations of vectorized array using the array+offset method
// with a few different lengths to test multiple variants. similar tests are
// vectorization_10 and vectorization_15

#include <deal.II/base/vectorization.h>

#include <limits>

#include "../tests.h"


template <typename Number, int n_numbers>
void
test()
{
  // since the number of array elements is system dependent, it is not a good
  // idea to print them to an output file. Instead, check the values manually
  const unsigned int      n_vectors = VectorizedArray<Number>::size();
  VectorizedArray<Number> arr[n_numbers];
  Number                  other[n_vectors * n_numbers];
  unsigned int            offsets[n_vectors];
  for (unsigned int v = 0; v < n_vectors; ++v)
    offsets[v] = v * n_numbers;

  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      other[i * n_numbers + j] = i * n_numbers + j + 13;

  vectorized_load_and_transpose(n_numbers, other, offsets, arr);
  unsigned int n_errors = 0;
  for (unsigned int j = 0; j < n_numbers; ++j)
    for (unsigned int i = 0; i < n_vectors; ++i)
      if (arr[j][i] != i * n_numbers + j + 13)
        ++n_errors;
  deallog << "load_and_transpose at          n=" << n_numbers
          << ": #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    for (unsigned int i = 0; i < n_numbers; ++i)
      {
        for (unsigned int j = 0; j < n_vectors; ++j)
          deallog << arr[i][j] << ' ';
        deallog << std::endl;
      }

  vectorized_transpose_and_store(true, n_numbers, arr, offsets, other);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      if (other[i * n_numbers + j] != 2. * (i * n_numbers + j + 13))
        ++n_errors;
  deallog << "transpose_and_store (  add) at n=" << n_numbers
          << ": #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    for (unsigned int i = 0; i < n_vectors; ++i)
      {
        for (unsigned int j = 0; j < n_numbers; ++j)
          deallog << other[i * n_numbers + j] << ' ';
        deallog << std::endl;
      }

  vectorized_transpose_and_store(false, n_numbers, arr, offsets, other);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      if (other[i * n_numbers + j] != (i * n_numbers + j + 13))
        ++n_errors;
  deallog << "transpose_and_store (noadd) at n=" << n_numbers
          << ": #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    for (unsigned int i = 0; i < n_vectors; ++i)
      {
        for (unsigned int j = 0; j < n_numbers; ++j)
          deallog << other[i * n_numbers + j] << ' ';
        deallog << std::endl;
      }
}



int
main()
{
  initlog();

  deallog.push("double");
  test<double, 1>();
  test<double, 2>();
  test<double, 3>();
  test<double, 4>();
  test<double, 5>();
  test<double, 6>();
  test<double, 7>();
  test<double, 8>();
  test<double, 9>();
  test<double, 10>();
  test<double, 32>();
  deallog.pop();
  deallog.push("float");
  test<float, 1>();
  test<float, 2>();
  test<float, 3>();
  test<float, 4>();
  test<float, 5>();
  test<float, 6>();
  test<float, 7>();
  test<float, 8>();
  test<float, 9>();
  test<float, 10>();
  test<float, 11>();
  test<float, 12>();
  test<float, 13>();
  test<float, 14>();
  test<float, 15>();
  test<float, 16>();
  test<float, 17>();
  test<float, 32>();
  deallog.pop();

  // test long double: in that case, the default
  // path of VectorizedArray is taken no matter
  // what was done for double or float
  deallog.push("long double");
  test<long double, 4>();
  deallog.pop();
}
