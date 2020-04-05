// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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
// for the set of all supported vectorization widths (otherwise the same as
// vectorization_05)

#include <deal.II/base/vectorization.h>

#include <limits>

#include "../tests.h"


template <typename Number, int n_numbers, int width>
void
do_test()
{
  // since the number of array elements is system dependent, it is not a good
  // idea to print them to an output file. Instead, check the values manually
  const unsigned int n_vectors = VectorizedArray<Number, width>::size();
  VectorizedArray<Number, width> arr[n_numbers];
  Number                         other[n_vectors * n_numbers];
  unsigned int                   offsets[n_vectors];
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
  if (n_errors > 0)
    {
      deallog << "load_and_transpose at          n=" << n_numbers
              << " width=" << width << ": #errors: " << n_errors << std::endl;

      for (unsigned int i = 0; i < n_numbers; ++i)
        {
          for (unsigned int j = 0; j < n_vectors; ++j)
            deallog << arr[i][j] << " ";
          deallog << std::endl;
        }
    }

  vectorized_transpose_and_store(true, n_numbers, arr, offsets, other);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      if (other[i * n_numbers + j] != 2. * (i * n_numbers + j))
        ++n_errors;
  if (n_errors > 0)
    {
      deallog << "transpose_and_store (  add) at n=" << n_numbers
              << " width=" << width << ": #errors: " << n_errors << std::endl;

      for (unsigned int i = 0; i < n_vectors; ++i)
        {
          for (unsigned int j = 0; j < n_numbers; ++j)
            deallog << other[i * n_numbers + j] << " ";
          deallog << std::endl;
        }
    }

  vectorized_transpose_and_store(false, n_numbers, arr, offsets, other);
  n_errors = 0;
  for (unsigned int i = 0; i < n_vectors; ++i)
    for (unsigned int j = 0; j < n_numbers; ++j)
      if (other[i * n_numbers + j] != (i * n_numbers + j))
        ++n_errors;
  if (n_errors > 0)
    {
      deallog << "transpose_and_store (noadd) at n=" << n_numbers
              << " width=" << width << ": #errors: " << n_errors << std::endl;

      for (unsigned int i = 0; i < n_vectors; ++i)
        {
          for (unsigned int j = 0; j < n_numbers; ++j)
            deallog << other[i * n_numbers + j] << " ";
          deallog << std::endl;
        }
    }
}


template <typename Number, int n_numbers>
struct Tester
{
  static void
  test()
  {
    do_test<Number, n_numbers, 1>();
  }
};

template <int n_numbers>
struct Tester<double, n_numbers>
{
  static void
  test()
  {
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
    do_test<double, n_numbers, 8>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
    do_test<double, n_numbers, 4>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
    do_test<double, n_numbers, 2>();
#endif

    do_test<double, n_numbers, 1>();
  }
};

template <int n_numbers>
struct Tester<float, n_numbers>
{
  static void
  test()
  {
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
    do_test<float, n_numbers, 16>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
    do_test<float, n_numbers, 8>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
    do_test<float, n_numbers, 4>();
#endif

    do_test<float, n_numbers, 1>();
  }
};


int
main()
{
  initlog();

  deallog.push("double");
  Tester<double, 1>::test();
  Tester<double, 9>::test();
  Tester<double, 32>::test();
  deallog.pop();
  deallog.push("float");
  Tester<float, 1>::test();
  Tester<float, 17>::test();
  Tester<float, 32>::test();
  deallog.pop();

  // test long double: in that case, the default
  // path of VectorizedArray is taken no matter
  // what was done for double or float
  deallog.push("long double");
  Tester<long double, 4>::test();
  deallog.pop();
}
