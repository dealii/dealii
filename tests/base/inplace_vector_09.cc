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


// Test inplace_vector's other free functions (erase() and erase_if())


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_erase()
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

  // erase()
  {
    std_cxx26::inplace_vector<T, 16> vec(as.begin(), as.end());
    for (unsigned int i = 0; i < 3; ++i)
      vec.insert(vec.end(), as.begin(), as.end());

    std_cxx26::erase(vec, as[1]);
    print(vec);
    deallog << std::endl;
  }

  // erase_if()
  {
    std_cxx26::inplace_vector<T, 16> vec(as.begin(), as.end());
    for (unsigned int i = 0; i < 3; ++i)
      vec.insert(vec.end(), as.begin(), as.end());

    std_cxx26::erase_if(vec, [&](const auto &a) {
      return a == as[0] || a == as[2];
    });
    print(vec);
    deallog << std::endl;
  }

  print_counts<T>();
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog << std::endl;
  deallog.push("A");
  test_erase<A>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_erase<int>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_erase<std::vector<int>>();
  deallog.pop();

  deallog.pop();
}
