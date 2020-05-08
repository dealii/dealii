// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// test printing a VectorizedArray object

#include <deal.II/base/vectorization.h>

#include <sstream>

#include "../tests.h"


template <typename Number>
void
test()
{
  constexpr unsigned int n = VectorizedArray<Number>::size();

  VectorizedArray<Number> a;
  std::stringstream       test_stream;
  for (unsigned int i = 0; i < n; ++i)
    a[i] = (i + 1);
  test_stream << a;

  std::stringstream reference_stream;
  for (unsigned int i = 0; i < n - 1; ++i)
    {
      reference_stream << Number(i + 1) << " ";
    }
  reference_stream << Number(n);

  if (reference_stream.str() != test_stream.str())
    {
      deallog << "Error: " << test_stream.str() << " should be "
              << reference_stream.str() << std::endl;
      AssertThrow(false, ExcInternalError());
    }
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  deallog << "Test double: ";
  test<double>();
  deallog << "Test int: ";
  test<int>();
}
