// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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
