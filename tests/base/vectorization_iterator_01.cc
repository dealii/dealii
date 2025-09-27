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


// test for arithmetic operations on VectorizedArray

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename Number>
void
test_const(const VectorizedArray<Number> &vector)
{
  unsigned int counter = 0;
  for (const auto v : vector)
    {
      AssertThrow(v == vector[counter++], ExcInternalError());
    }
}

template <typename Number>
void
test_nonconst(VectorizedArray<Number> &vector)
{
  unsigned int counter = 0;
  for (auto v : vector)
    {
      AssertThrow(v == vector[counter++], ExcInternalError());
    }
}

template <typename Number>
void
test()
{
  VectorizedArray<Number> vector;

  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    vector[v] = v;

  test_const(vector);
  test_nonconst(vector);
}


int
main()
{
  initlog();

  deallog.push("float");
  test<float>();
  deallog << "OK" << std::endl;
  deallog.pop();

  deallog.push("double");
  test<double>();
  deallog << "OK" << std::endl;
  deallog.pop();
}
