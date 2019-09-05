// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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


// test Utilities::fixed_power for vectorized array

#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"

template <typename VectorizedArrayType>
void
do_test(const VectorizedArrayType array)
{
  deallog << "  test " << VectorizedArrayType::n_array_elements
          << " array elements" << std::endl;

  auto exponentiated_array = Utilities::fixed_power<3>(array);

  for (unsigned int i = 0; i < VectorizedArrayType::n_array_elements; i++)
    deallog << exponentiated_array[i] << " ";
  deallog << std::endl;

  exponentiated_array = Utilities::fixed_power<-3>(array);

  for (unsigned int i = 0; i < VectorizedArrayType::n_array_elements; i++)
    deallog << exponentiated_array[i] << " ";
  deallog << std::endl;
}


int
main()
{
  initlog();

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 3 && defined(__AVX512F__)
  do_test(VectorizedArray<double, 8>(2.0));
  do_test(VectorizedArray<float, 16>(2.0));
#endif

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 && defined(__AVX__)
  do_test(VectorizedArray<double, 4>(2.0));
  do_test(VectorizedArray<float, 8>(2.0));
#endif

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 1 && defined(__SSE2__)
  do_test(VectorizedArray<double, 2>(2.0));
  do_test(VectorizedArray<float, 4>(2.0));
#endif

  do_test(VectorizedArray<double, 1>(2.0));
  do_test(VectorizedArray<float, 1>(2.0));
}
