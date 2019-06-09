// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// tests column padding for different sizes

#include <deal.II/lac/block_csr_matrix.h>

#include <fstream>
#include <iostream>

using namespace dealii;



template <typename NumberType = double>
void
test(const std::vector<unsigned int> sizes)
{
  const NumberType   dummy = NumberType(0.);
  const unsigned int CL    = 64;
  std::cout << "Testing  " << dealii::Utilities::type_to_string(dummy)
            << std::endl
            << "  size:  " << sizeof(NumberType) << std::endl
            << "  SIMD:  " << VectorizedArray<NumberType>::n_array_elements
            << std::endl
            << "  CL:    " << CL << std::endl
            << "  el/CL: " << CL / sizeof(NumberType) << std::endl
            << std::endl;

  for (const auto s : sizes)
    {
      const auto upper = internal::padded_size<NumberType>(s);
      Assert(internal::ceil_divisible_by(
               s, VectorizedArray<NumberType>::n_array_elements) <= upper,
             ExcInternalError());
      Assert(internal::ceil_divisible_by(s, CL / sizeof(NumberType)) == upper,
             ExcInternalError());
      std::cout << s << " -> " << upper << std::endl;
    }
}

int
main()
{
  test<double>({3, 4, 7, 8, 12, 15, 16, 24, 32});
  test<float>({3, 4, 7, 8, 12, 15, 16, 24, 32});

  return 0;
}
