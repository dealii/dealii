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


// Test std_cxx26::inplace_vector's constructors.


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_ctors()
{
  constexpr bool is_a = std::is_same_v<T, A>;

  deallog << "default ctor" << std::endl;
  {
    std_cxx26::inplace_vector<T, 16> vec;
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "count default ctor" << std::endl;
  {
    std_cxx26::inplace_vector<T, 16> vec(2);
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "count copy ctor" << std::endl;
  {
    T a{};
    if constexpr (std::is_same_v<T, int>)
      a = 42;
    if constexpr (std::is_same_v<T, std::vector<int>>)
      a = {3, 5};
    std_cxx26::inplace_vector<T, 16> vec(2, a);
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
  deallog << std::endl;

  deallog << "iterator range, copy, and move ctor" << std::endl;
  {
    std::array<T, 2> as{};
    if constexpr (std::is_same_v<T, int>)
      std::iota(as.begin(), as.end(), -11);
    if constexpr (std::is_same_v<T, std::vector<int>>)
      {
        as[0] = {1, 1, 2};
        as[1] = {3, 5};
      }
    std_cxx26::inplace_vector<T, 16> vec1(as.begin(), as.end());
    print(vec1);
    deallog << std::endl;
    deallog << "copy ctor" << std::endl;
    std_cxx26::inplace_vector<T, 16> vec2(vec1);
    print(vec2);
    deallog << std::endl;
    deallog << "and move ctor" << std::endl;
    std_cxx26::inplace_vector<T, 16> vec3(std::move(vec1));
    deallog << "moved-from:" << std::endl;
    print(vec1);
    deallog << std::endl;
    deallog << "moved-to:" << std::endl;
    print(vec3);
    deallog << std::endl;
  }
  print_counts<T>();

  deallog << "initializer_list ctor" << std::endl;
  {
    std::array<T, 2> as{};
    if constexpr (std::is_same_v<T, int>)
      std::iota(as.begin(), as.end(), -11);
    if constexpr (std::is_same_v<T, std::vector<int>>)
      {
        as[0] = {1, 1, 2};
        as[1] = {3, 5};
      }
    std_cxx26::inplace_vector<T, 16> vec(
      std::initializer_list<T>{as[0], as[1]});
    print(vec);
    deallog << std::endl;
  }
  print_counts<T>();
}

int
main()
{
  initlog();

  deallog.push("ctors");
  deallog << std::endl;
  deallog.push("A");
  test_ctors<A>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_ctors<int>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_ctors<std::vector<int>>();
  deallog.pop();

  deallog.pop();
}
