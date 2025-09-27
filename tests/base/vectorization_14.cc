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


// test compare_and_apply_mask

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename VectorizedArrayType>
void
do_test()
{
  deallog << "  test " << VectorizedArrayType::size() << " array elements"
          << std::endl;

  VectorizedArrayType left;
  for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
    left[i] = i + 1.;

  VectorizedArrayType right(3.);
  VectorizedArrayType true_values(1.);
  VectorizedArrayType false_values(-1.);

  auto result = compare_and_apply_mask<SIMDComparison::equal>(left,
                                                              right,
                                                              true_values,
                                                              false_values);
  deallog << result << std::endl;

  result = compare_and_apply_mask<SIMDComparison::not_equal>(left,
                                                             right,
                                                             true_values,
                                                             false_values);
  deallog << result << std::endl;

  result = compare_and_apply_mask<SIMDComparison::less_than>(left,
                                                             right,
                                                             true_values,
                                                             false_values);
  deallog << result << std::endl;

  result = compare_and_apply_mask<SIMDComparison::less_than_or_equal>(
    left, right, true_values, false_values);
  deallog << result << std::endl;

  result = compare_and_apply_mask<SIMDComparison::greater_than>(left,
                                                                right,
                                                                true_values,
                                                                false_values);
  deallog << result << std::endl;

  result = compare_and_apply_mask<SIMDComparison::greater_than_or_equal>(
    left, right, true_values, false_values);
  deallog << result << std::endl;
}


int
main()
{
  initlog();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test<VectorizedArray<double, 8>>();
  do_test<VectorizedArray<float, 16>>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test<VectorizedArray<double, 4>>();
  do_test<VectorizedArray<float, 8>>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test<VectorizedArray<double, 2>>();
  do_test<VectorizedArray<float, 4>>();
#endif

  do_test<VectorizedArray<double, 1>>();
  do_test<VectorizedArray<float, 1>>();
}
