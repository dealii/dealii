// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test constructor of VectorizedArray with the given initializer list.

#include <deal.II/base/vectorization.h>

#include <limits>

#include "../tests.h"


template <typename Number, int width>
void
do_test(const std::initializer_list<Number> &list)
{
  const unsigned int entries = std::distance(list.begin(), list.end());

  if (entries > width)
    return;

  const VectorizedArray<Number, width> vec = list;


  for (unsigned int i = 0; i < entries; ++i)
    AssertDimension(vec[i], i + 1);

  for (unsigned int i = entries; i < width; ++i)
    AssertDimension(vec[i], 0);
}

template <typename Number, int width>
void
do_test()
{
  do_test<Number, width>({1});
  do_test<Number, width>({1, 2});
  do_test<Number, width>({1, 2, 3});
  do_test<Number, width>({1, 2, 3, 4});
  do_test<Number, width>({1, 2, 3, 4, 5});
  do_test<Number, width>({1, 2, 3, 4, 5, 6});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
  do_test<Number, width>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
  do_test<Number, width>(
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16});
}


int
main()
{
  initlog();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test<float, 16>();
  do_test<double, 8>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test<float, 8>();
  do_test<double, 4>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test<float, 4>();
  do_test<double, 2>();
#endif

  do_test<float, 1>();
  do_test<double, 1>();

  deallog << "OK!" << std::endl;
}
