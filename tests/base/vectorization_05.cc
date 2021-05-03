// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


// test transpose operations of vectorized array using the array+offset method
// (otherwise the same as vectorization_05)

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
      other[i * n_numbers + j] = i * n_numbers + j;

  vectorized_load_and_transpose(n_numbers, other, offsets, arr);
  unsigned int n_errors = 0;
  for (unsigned int j = 0; j < n_numbers; ++j)
    for (unsigned int i = 0; i < n_vectors; ++i)
      if (arr[j][i] != i * n_numbers + j)
        ++n_errors;
  deallog << "load_and_transpose at          n=" << n_numbers
          << ": #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    for (unsigned int i = 0; i < n_numbers; ++i)
      {
        for (unsigned int j = 0; j < n_vectors; ++j)
          deallog << arr[i][j] << " ";
        deallog << std::endl;
      }

  vectorized_transpose_and_store(true, n_numbers, arr, offsets, other);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      if (other[i * n_numbers + j] != 2. * (i * n_numbers + j))
        ++n_errors;
  deallog << "transpose_and_store (  add) at n=" << n_numbers
          << ": #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    for (unsigned int i = 0; i < n_vectors; ++i)
      {
        for (unsigned int j = 0; j < n_numbers; ++j)
          deallog << other[i * n_numbers + j] << " ";
        deallog << std::endl;
      }

  vectorized_transpose_and_store(false, n_numbers, arr, offsets, other);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      if (other[i * n_numbers + j] != (i * n_numbers + j))
        ++n_errors;
  deallog << "transpose_and_store (noadd) at n=" << n_numbers
          << ": #errors: " << n_errors << std::endl;
  if (n_errors > 0)
    for (unsigned int i = 0; i < n_vectors; ++i)
      {
        for (unsigned int j = 0; j < n_numbers; ++j)
          deallog << other[i * n_numbers + j] << " ";
        deallog << std::endl;
      }
}



int
main()
{
  initlog();

  deallog.push("double");
  test<double, 1>();
  test<double, 9>();
  test<double, 32>();
  deallog.pop();
  deallog.push("float");
  test<float, 1>();
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
