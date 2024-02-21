// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test std::pow for vectorized array

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename VectorizedArrayType>
void
do_test(const VectorizedArrayType                      array,
        const typename VectorizedArrayType::value_type number)
{
  deallog << "  test " << VectorizedArrayType::size() << " array elements"
          << std::endl;

  auto exponentiated_array = std::pow(array, number);

  for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
    deallog << exponentiated_array[i] << ' ';
  deallog << std::endl;
}


template <typename VectorizedArrayType>
void
do_test(const VectorizedArrayType array, const VectorizedArrayType number)
{
  deallog << "  test " << VectorizedArrayType::size() << " array elements"
          << std::endl;

  auto exponentiated_array = std::pow(array, number);

  for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
    deallog << exponentiated_array[i] << ' ';
  deallog << std::endl;
}


int
main()
{
  initlog();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test(VectorizedArray<double, 8>(2.0), 3.0);
  do_test(VectorizedArray<float, 16>(2.0), 3.0);

  do_test(VectorizedArray<double, 8>(2.0), VectorizedArray<double, 8>(3.0));
  do_test(VectorizedArray<float, 16>(2.0), VectorizedArray<float, 16>(3.0));
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test(VectorizedArray<double, 4>(2.0), 3.0);
  do_test(VectorizedArray<float, 8>(2.0), 3.0);

  do_test(VectorizedArray<double, 4>(2.0), VectorizedArray<double, 4>(3.0));
  do_test(VectorizedArray<float, 8>(2.0), VectorizedArray<float, 8>(3.0));
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test(VectorizedArray<double, 2>(2.0), 3.0);
  do_test(VectorizedArray<float, 4>(2.0), 3.0);

  do_test(VectorizedArray<double, 2>(2.0), VectorizedArray<double, 2>(3.0));
  do_test(VectorizedArray<float, 4>(2.0), VectorizedArray<float, 4>(3.0));
#endif

  do_test(VectorizedArray<double, 1>(2.0), 3.0);
  do_test(VectorizedArray<float, 1>(2.0), 3.0);

  do_test(VectorizedArray<double, 1>(2.0), VectorizedArray<double, 1>(3.0));
  do_test(VectorizedArray<float, 1>(2.0), VectorizedArray<float, 1>(3.0));
}
