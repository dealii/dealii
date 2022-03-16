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


// test mixed precision std::pow for vectorized array

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename VectorizedArrayBaseType, typename ExponentType>
void
do_test(const VectorizedArrayBaseType array, const ExponentType number)
{
  deallog << "  test " << VectorizedArrayBaseType::size() << " array elements"
          << std::endl;

  const auto exponentiated_array = std::pow(array, number);

  for (unsigned int i = 0; i < VectorizedArrayBaseType::size(); ++i)
    deallog << exponentiated_array[i] << ' ';
  deallog << std::endl;
}


int
main()
{
  initlog();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test(VectorizedArray<double, 8>(2.0), 3.0f);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test(VectorizedArray<double, 4>(2.0), 3.0f);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test(VectorizedArray<double, 2>(2.0), 3.0f);
#endif

  do_test(VectorizedArray<double, 1>(2.0), 3.0);
  do_test(VectorizedArray<float, 1>(2.0f), 3.0);

  do_test(VectorizedArray<double, 1>(2.0), VectorizedArray<float, 1>(3.0f));
  do_test(VectorizedArray<float, 1>(2.0f), VectorizedArray<double, 1>(3.0));
}
