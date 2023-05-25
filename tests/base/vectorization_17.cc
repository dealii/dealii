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


// test mixed precision arithmetic for vectorized array

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename VectorizedArrayType>
void
print(const std::string &msg, const VectorizedArrayType &array)
{
  deallog << msg << "  ";
  for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
    deallog << array[i] << ' ';
  deallog << std::endl;
}


template <typename VectorizedArrayType, typename NumberType2>
void
do_test(const VectorizedArrayType array, const NumberType2 number)
{
  deallog << "  test " << VectorizedArrayType::size() << " array elements"
          << std::endl;

  const auto addition_1 = array + number;
  const auto addition_2 = number + array;

  print("add 1", addition_1);
  print("add 2", addition_2);

  const auto subtraction_1 = array - number;
  const auto subtraction_2 = number - array;

  print("sub 1", subtraction_1);
  print("sub 2", subtraction_2);

  const auto multiplication_1 = array * number;
  const auto multiplication_2 = number * array;

  print("mult 1", multiplication_1);
  print("mult 2", multiplication_2);

  const auto division_1 = array / number;
  const auto division_2 = number / array;

  print("div 1", division_1);
  print("div 2", division_2);
}


int
main()
{
  initlog();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test(VectorizedArray<double, 8>(2.0), 3.0f);
  do_test(VectorizedArray<float, 16>(2.0), 3.0);

  do_test(VectorizedArray<double, 8>(2.0), VectorizedArray<float, 16>(3.0));
  do_test(VectorizedArray<float, 16>(2.0), VectorizedArray<double, 8>(3.0));
  deallog << std::endl;
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test(VectorizedArray<double, 4>(2.0), 3.0f);
  do_test(VectorizedArray<float, 8>(2.0), 3.0);

  do_test(VectorizedArray<double, 4>(2.0), VectorizedArray<float, 8>(3.0));
  do_test(VectorizedArray<float, 8>(2.0), VectorizedArray<double, 4>(3.0));
  deallog << std::endl;
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test(VectorizedArray<double, 2>(2.0), 3.0f);
  do_test(VectorizedArray<float, 4>(2.0), 3.0);

  do_test(VectorizedArray<double, 2>(2.0), VectorizedArray<float, 4>(3.0));
  do_test(VectorizedArray<float, 4>(2.0), VectorizedArray<double, 2>(3.0));
  deallog << std::endl;
#endif

  do_test(VectorizedArray<double, 1>(2.0), 3.0f);
  do_test(VectorizedArray<double, 1>(2.0), std::complex<double>{3.0, 0.0});

  do_test(VectorizedArray<float, 1>(2.0), 3.0f);
  do_test(VectorizedArray<float, 1>(2.0), std::complex<float>{3.0f, 0.0f});

  do_test(VectorizedArray<double, 1>(2.0), 3.0f);
  do_test(VectorizedArray<double, 1>(2.0), std::complex<float>{3.0f, 0.0f});
}
