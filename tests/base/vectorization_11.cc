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


// test make_vectorized_array function for all possible value_types and vector
// lengths

#include <deal.II/base/vectorization.h>

#include <limits>

#include "../tests.h"


template <typename VectorizedArrayType>
void
do_test(const VectorizedArrayType                      array,
        const typename VectorizedArrayType::value_type number)
{
  deallog << "  test " << VectorizedArrayType::size() << " array elements"
          << std::endl;
  for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
    if (array[i] != number)
      deallog << "  problem in element " << i << std::endl;
}


template <typename Number>
struct Tester
{
  static void
  test()
  {}
};

template <>
struct Tester<double>
{
  static void
  test()
  {
    do_test(make_vectorized_array<double>(2.0), 2.0);
    do_test(VectorizedArray<double>(2.0), 2.0);

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
    do_test(make_vectorized_array<VectorizedArray<double, 8>>(2.0), 2.0);
    do_test(VectorizedArray<double, 8>(2.0), 2.0);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
    do_test(make_vectorized_array<VectorizedArray<double, 4>>(2.0), 2.0);
    do_test(VectorizedArray<double, 4>(2.0), 2.0);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
    do_test(make_vectorized_array<VectorizedArray<double, 2>>(2.0), 2.0);
    do_test(VectorizedArray<double, 2>(2.0), 2.0);
#endif

    do_test(make_vectorized_array<VectorizedArray<double, 1>>(2.0), 2.0);
    do_test(VectorizedArray<double, 1>(2.0), 2.0);
  }
};

template <>
struct Tester<float>
{
  static void
  test()
  {
    do_test(make_vectorized_array<float>(2.0), 2.0);
    do_test(VectorizedArray<float>(2.0), 2.0);

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
    do_test(make_vectorized_array<VectorizedArray<float, 16>>(2.0), 2.0);
    do_test(VectorizedArray<float, 16>(2.0), 2.0);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
    do_test(make_vectorized_array<VectorizedArray<float, 8>>(2.0), 2.0);
    do_test(VectorizedArray<float, 8>(2.0), 2.0);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
    do_test(make_vectorized_array<VectorizedArray<float, 4>>(2.0), 2.0);
    do_test(VectorizedArray<float, 4>(2.0), 2.0);
#endif

    do_test(make_vectorized_array<VectorizedArray<float, 1>>(2.0), 2.0);
    do_test(VectorizedArray<float, 1>(2.0), 2.0);
  }
};


int
main()
{
  initlog();

  deallog << "double:" << std::endl;
  Tester<double>::test();

  deallog << "float:" << std::endl;
  Tester<float>::test();
}
