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


// Test std_cxx26::inplace_vector's insertion functions (emplace(), insert(),
// erase(), clear(), and swap())


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_access()
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

  // emplace
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.emplace(vec.cbegin() + 1);
    print(vec);
    deallog << std::endl;
    vec.emplace(vec.cbegin() + 1, as[2]);
    print(vec);
    deallog << std::endl;
    if constexpr (std::is_same_v<T, std::vector<int>>)
      {
        vec.emplace(vec.cbegin() + 1, 5u, 42);
        print(vec);
        deallog << std::endl;
      }
  }

  // insert value, insert move, and insert count and value
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.insert(vec.cbegin() + 1, as[0]);
    print(vec);
    deallog << std::endl;

    auto copy = as[1];
    vec.insert(vec.cbegin() + 1, std::move(copy));
    if constexpr (std::is_same_v<T, std::vector<int>>)
      AssertThrow(copy.empty(), ExcInternalError());
    print(vec);
    deallog << std::endl;

    auto copy2 = as[2];
    vec.insert(vec.cend(), 3, copy2);
    if constexpr (std::is_same_v<T, std::vector<int>>)
      AssertThrow(copy2 == as[2], ExcInternalError());
    print(vec);
    deallog << std::endl;
  }

  // insert range
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.insert(vec.cbegin() + 1, as.rbegin(), as.rend());
    print(vec);
    deallog << std::endl;
  }

  // insert initializer list
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.insert(vec.cbegin() + 1, {as[0], as[2], as[1]});
    print(vec);
    deallog << std::endl;
  }

  // erase position
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.erase(vec.cend() - 1);
    vec.erase(vec.cbegin());
    print(vec);
    deallog << std::endl;
  }

  // erase iterator range
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.erase(vec.begin() + 1, vec.end());
    print(vec);
    deallog << std::endl;
  }

  // clear
  {
    std_cxx26::inplace_vector<T, 8> vec(as.begin(), as.end());

    vec.clear();
    print(vec);
    deallog << std::endl;
  }

  // swap
  {
    std_cxx26::inplace_vector<T, 8> vec1(as.begin(), as.end()), vec2(5);
    vec1.insert(vec1.end(), as.begin(), as.end());

    vec1.swap(vec2);
    print(vec1);
    deallog << std::endl;

    using std::swap;
    swap(vec1, vec2);
    print(vec1);
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
  test_access<A>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_access<int>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_access<std::vector<int>>();
  deallog.pop();

  deallog.pop();
}
