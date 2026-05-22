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


// Test inplace_vector's comparison operators.


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T>
void
test_compare(const bool different_sizes)
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

  std_cxx26::inplace_vector<T, 16> vec(as.begin(), as.end()), vec2, vec3;
  for (unsigned int i = 0; i < 3; ++i)
    vec.insert(vec.end(), as.begin(), as.end());
  vec2 = vec;
  vec3 = vec;
  if (different_sizes)
    {
      vec2.resize(vec2.size() - 4);
      vec3.push_back(as[0]);
    }

  std::sort(vec2.begin(), vec2.end());
  std::sort(vec3.rbegin(), vec3.rend());

  print(vec);
  deallog << std::endl;
  print(vec2);
  deallog << std::endl;
  print(vec3);
  deallog << std::endl;

  deallog << "vec == vec  = " << (vec == vec) << std::endl;
  deallog << "vec != vec  = " << (vec != vec) << std::endl;
  deallog << "vec < vec   = " << (vec < vec) << std::endl;
  deallog << "vec > vec   = " << (vec > vec) << std::endl;
  deallog << "vec <= vec  = " << (vec <= vec) << std::endl;
  deallog << "vec >= vec  = " << (vec >= vec) << std::endl;

  deallog << "vec == vec2 = " << (vec == vec2) << std::endl;
  deallog << "vec != vec2 = " << (vec != vec2) << std::endl;
  deallog << "vec < vec2  = " << (vec < vec2) << std::endl;
  deallog << "vec > vec2  = " << (vec > vec2) << std::endl;
  deallog << "vec <= vec2 = " << (vec <= vec2) << std::endl;
  deallog << "vec >= vec2 = " << (vec >= vec2) << std::endl;

  deallog << "vec == vec3 = " << (vec == vec3) << std::endl;
  deallog << "vec != vec3 = " << (vec != vec3) << std::endl;
  deallog << "vec < vec3  = " << (vec < vec3) << std::endl;
  deallog << "vec > vec3  = " << (vec > vec3) << std::endl;
  deallog << "vec <= vec3 = " << (vec <= vec3) << std::endl;
  deallog << "vec >= vec3 = " << (vec >= vec3) << std::endl;

  deallog << std::endl;
}

int
main()
{
  initlog();

  A::logging() = false;
  deallog << std::endl;
  deallog.push("A");
  test_compare<A>(false);
  test_compare<A>(true);
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_compare<int>(false);
  test_compare<int>(true);
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_compare<std::vector<int>>(false);
  test_compare<std::vector<int>>(true);
  deallog.pop();

  deallog.pop();
}
