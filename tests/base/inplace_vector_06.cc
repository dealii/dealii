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


// Test std_cxx26::inplace_vector's capacity (empty(), size(), max_size(), etc.)
// functions.


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_capacity()
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

  std_cxx26::inplace_vector<T, 4> vec;

  auto print_sizes = [&]() {
    deallog << "empty    = " << vec.empty() << std::endl;
    deallog << "size     = " << vec.size() << std::endl;
    deallog << "max_size = " << vec.max_size() << std::endl;
    deallog << "capacity = " << vec.capacity() << std::endl;
  };

  print_sizes();
  print(vec);
  deallog << std::endl;

  vec.resize(2);
  print_sizes();
  print(vec);
  deallog << std::endl;

  vec.resize(3, as[0]);
  print_sizes();
  print(vec);
  deallog << std::endl;

  // test no-op functions
  const auto vec2 = vec;
  vec2.shrink_to_fit();
  AssertThrow(vec2 == vec, ExcInternalError());
  vec2.reserve(4);
  AssertThrow(vec2 == vec, ExcInternalError());

  print_counts<T>();
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog << std::endl;
  deallog.push("A");
  test_capacity<A>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_capacity<int>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_capacity<std::vector<int>>();
  deallog.pop();

  deallog.pop();
}
