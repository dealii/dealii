// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for std::max_element and std::distance on VectorizedArray

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename Number>
void
test_const(const VectorizedArray<Number> &vector)
{
  AssertDimension(*std::max_element(vector.begin(), vector.end()),
                  VectorizedArray<Number>::size() - 1);
  AssertDimension(std::distance(vector.begin(),
                                std::max_element(vector.begin(), vector.end())),
                  VectorizedArray<Number>::size() - 1);
  AssertDimension(std::distance(vector.begin(),
                                vector.begin() +
                                  (VectorizedArray<Number>::size() - 1)),
                  VectorizedArray<Number>::size() - 1);
}

template <typename Number>
void
test_nonconst(VectorizedArray<Number> &vector)
{
  AssertDimension(*std::max_element(vector.begin(), vector.end()),
                  VectorizedArray<Number>::size() - 1);
  AssertDimension(std::distance(vector.begin(),
                                std::max_element(vector.begin(), vector.end())),
                  VectorizedArray<Number>::size() - 1);
  AssertDimension(std::distance(vector.begin(),
                                vector.begin() +
                                  (VectorizedArray<Number>::size() - 1)),
                  VectorizedArray<Number>::size() - 1);

  auto it = vector.begin();
  std::advance(it, VectorizedArray<Number>::size() - 1);
  AssertDimension(*it, VectorizedArray<Number>::size() - 1);
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
