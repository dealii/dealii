// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test std_cxx26::inplace_vector's iterator (begin(), end(), etc.) interface.


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_iterators()
{
  constexpr bool is_a = std::is_same_v<T, A>;

  std::array<T, 3> as{};
  if constexpr (std::is_same_v<T, int>)
    std::iota(as.begin(), as.end(), -11);
  if constexpr (std::is_same_v<T, std::vector<int>>)
    {
      as[0] = {1, 1, 2};
      as[1] = {42};
      as[2] = {3, 5};
    }

  std_cxx26::inplace_vector<T, 16> vec;
  const auto                      &cvec = vec;
  vec.assign(as.begin(), as.end());
  print(vec);
  deallog << std::endl;
  AssertThrow(vec.begin() + vec.size() == vec.end(), ExcInternalError());
  AssertThrow(cvec.begin() + cvec.size() == cvec.end(), ExcInternalError());
  AssertThrow(vec.cbegin() + vec.size() == vec.cend(), ExcInternalError());
  AssertThrow(vec.rbegin() + vec.size() == vec.rend(), ExcInternalError());
  AssertThrow(cvec.rbegin() + vec.size() == cvec.rend(), ExcInternalError());
  AssertThrow(vec.crbegin() + vec.size() == vec.crend(), ExcInternalError());

  std::size_t i = 0;
  for (auto it = vec.begin(); it != vec.end(); ++it)
    {
      AssertThrow(*it == as[i], ExcInternalError());
      AssertThrow(std::addressof(*it) == std::addressof(vec[i]),
                  ExcInternalError());
      ++i;
    }

  i = 0;
  for (auto it = cvec.begin(); it != cvec.end(); ++it)
    {
      AssertThrow(*it == as[i], ExcInternalError());
      AssertThrow(std::addressof(*it) == std::addressof(cvec[i]),
                  ExcInternalError());
      ++i;
    }

  i = vec.size() - 1;
  for (auto it = vec.rbegin(); it != vec.rend(); ++it)
    {
      AssertThrow(*it == as[i], ExcInternalError());
      AssertThrow(std::addressof(*it) == std::addressof(vec[i]),
                  ExcInternalError());
      --i;
    }

  i = vec.size() - 1;
  for (auto it = cvec.rbegin(); it != cvec.rend(); ++it)
    {
      AssertThrow(*it == as[i], ExcInternalError());
      AssertThrow(std::addressof(*it) == std::addressof(cvec[i]),
                  ExcInternalError());
      --i;
    }

  i = vec.size() - 1;
  for (auto it = vec.crbegin(); it != vec.crend(); ++it)
    {
      AssertThrow(*it == as[i], ExcInternalError());
      AssertThrow(std::addressof(*it) == std::addressof(cvec[i]),
                  ExcInternalError());
      --i;
    }

  std::sort(vec.begin(), vec.end());
  print(vec);
  deallog << std::endl;

  std::sort(vec.rbegin(), vec.rend());
  print(vec);
  deallog << std::endl;

  print_counts<T>();
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog << std::endl;
  deallog.push("A");
  test_iterators<A>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_iterators<int>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_iterators<std::vector<int>>();
  deallog.pop();

  deallog.pop();
}
