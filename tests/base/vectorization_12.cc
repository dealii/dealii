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
      reference_stream << Number(i + 1) << ' ';
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
