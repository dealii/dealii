// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


// test AlignedVector::data

#include <deal.II/base/aligned_vector.h>

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
  AssertThrow(array.data() == &array[0], ExcInternalError());
}



template <typename T>
void
test()
{
  AlignedVector<T>        array;
  const AlignedVector<T> &const_array = array;
  test_array_empty(array);
  test_array_empty(const_array);

  array.push_back(T{});
  test_array_single(array);
  test_array_single(const_array);
}



int
main()
{
  initlog();

  test<int>();
  test<long double>();
  test<AlignedVector<double>>();

  deallog << "OK" << std::endl;
}
