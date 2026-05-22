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


// Test std_cxx26::inplace_vector's emplace_back(), push_back(), and pop_back()
// functions. We can't use try_emplace_back etc. since those return
// std::optional<T&>, which was not supported before C++26.


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_emplace_pop_push()
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

  deallog << "emplace_back" << std::endl;
  {
    std_cxx26::inplace_vector<T, 16> vec;
    vec.assign(as.begin(), as.end());
    print(vec);
    deallog << std::endl;
    auto copy_as = as;
    vec.emplace_back(as[0]);
    vec.emplace_back(std::move(copy_as[1]));
    if constexpr (std::is_same_v<T, std::vector<int>>)
      {
        vec.emplace_back(5u, 42);
        AssertThrow(copy_as[0].size() == 3, ExcInternalError());
        AssertThrow(copy_as[1].size() == 0, ExcInternalError());
      }
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "push_back" << std::endl;
  {
    std_cxx26::inplace_vector<T, 16> vec;
    vec.assign(as.begin(), as.end());
    print(vec);
    deallog << std::endl;
    auto copy_as = as;
    vec.push_back(copy_as[0]);
    vec.push_back(std::move(copy_as[1]));
    if constexpr (std::is_same_v<T, std::vector<int>>)
      {
        AssertThrow(copy_as[0].size() == 3, ExcInternalError());
        AssertThrow(copy_as[1].size() == 0, ExcInternalError());
      }
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "pop_back" << std::endl;
  {
    std_cxx26::inplace_vector<T, 16> vec;
    vec.assign(as.begin(), as.end());
    print(vec);
    deallog << std::endl;
    vec.pop_back();
    AssertThrow(vec.size() == 2, ExcInternalError());
    print(vec);
    deallog << std::endl;
    vec.pop_back();
    print(vec);
    AssertThrow(vec.size() == 1, ExcInternalError());
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;


  deallog << "unchecked_emplace_back" << std::endl;
  {
    std_cxx26::inplace_vector<T, 16> vec;
    vec.assign(as.begin(), as.end());
    print(vec);
    deallog << std::endl;
    auto copy_as = as;
    vec.unchecked_emplace_back(copy_as[0]);
    vec.unchecked_emplace_back(std::move(copy_as[1]));
    if constexpr (std::is_same_v<T, std::vector<int>>)
      {
        vec.unchecked_emplace_back(5u, 42);
        AssertThrow(copy_as[0].size() == 3, ExcInternalError());
        AssertThrow(copy_as[1].size() == 0, ExcInternalError());
      }
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "unchecked_push_back (move)" << std::endl;
  {
    std_cxx26::inplace_vector<T, 6> vec;
    vec.assign(as.begin(), as.end());
    print(vec);
    deallog << std::endl;
    auto copy_as = as;

    for (std::size_t i = 0; i < 3; ++i)
      {
        auto &back =
          vec.unchecked_push_back(std::move(copy_as[i % copy_as.size()]));
        if constexpr (std::is_same_v<T, std::vector<int>>)
          AssertThrow(copy_as[i % copy_as.size()].size() == 0,
                      ExcInternalError());
        AssertThrow(&back == &vec.back(), ExcInternalError());
        print(vec);
        deallog << std::endl;
      }
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "unchecked_push_back (copy)" << std::endl;
  {
    std_cxx26::inplace_vector<T, 6> vec;
    vec.assign(as.begin(), as.end());
    print(vec);
    deallog << std::endl;
    auto copy_as = as;

    for (std::size_t i = 0; i < 3; ++i)
      {
        auto &back = vec.unchecked_push_back(copy_as[i % copy_as.size()]);
        AssertThrow(&back == &vec.back(), ExcInternalError());
        print(vec);
        deallog << std::endl;
      }
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
  test_emplace_pop_push<A>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_emplace_pop_push<int>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_emplace_pop_push<std::vector<int>>();
  deallog.pop();

  deallog.pop();
}
