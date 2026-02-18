// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test Vector::data

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector.templates.h>

#include "../tests.h"


template <typename Container>
void
test_array_empty(Container &array)
{
  AssertThrow(array.begin() == array.end(), ExcInternalError());
  AssertThrow(array.end() == array.data(), ExcInternalError());
  AssertThrow(array.data() == nullptr, ExcInternalError());
}



template <typename Container>
void
test_array_single(Container &array)
{
  AssertThrow(array.begin() + 1 == array.end(), ExcInternalError());
  AssertThrow(array.end() - 1 == array.data(), ExcInternalError());
}



template <typename T>
void
test()
{
  Vector<T>        array;
  const Vector<T> &const_array = array;
  test_array_empty(array);
  test_array_empty(const_array);

  array.reinit(1);
  test_array_single(array);
  test_array_single(const_array);
}



int
main()
{
  initlog();

  test<int>();
  test<double>();
  test<long double>();

  deallog << "OK" << std::endl;
}
