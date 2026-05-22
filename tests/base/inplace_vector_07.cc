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


// Test std_cxx26::inplace_vector's access operators (operator[], at(), front(),
// etc.)


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

  std_cxx26::inplace_vector<T, 4> vec(as.begin(), as.end());
  const auto                     &cvec = vec;

  for (std::size_t i = 0; i < as.size(); ++i)
    {
      deallog << "vec[" << i << "] = " << vec[i] << std::endl;
      AssertThrow(std::addressof(vec[i]) == std::addressof(vec.at(i)),
                  ExcInternalError());
    }
  AssertThrow(std::addressof(vec[0]) == std::addressof(vec.front()),
              ExcInternalError());
  AssertThrow(std::addressof(vec[2]) == std::addressof(vec.back()),
              ExcInternalError());
  AssertThrow(std::addressof(vec[0]) == vec.data(), ExcInternalError());

  try
    {
      vec.at(42);
    }
  catch (const std::out_of_range &)
    {
      deallog << "caught correct exception" << std::endl;
    }

  deallog << "const versions:" << std::endl;
  for (std::size_t i = 0; i < as.size(); ++i)
    {
      deallog << "vec[" << i << "] = " << cvec[i] << std::endl;
      AssertThrow(std::addressof(cvec[i]) == std::addressof(cvec.at(i)),
                  ExcInternalError());
    }
  AssertThrow(std::addressof(cvec[0]) == std::addressof(cvec.front()),
              ExcInternalError());
  AssertThrow(std::addressof(cvec[2]) == std::addressof(cvec.back()),
              ExcInternalError());
  AssertThrow(std::addressof(cvec[0]) == cvec.data(), ExcInternalError());

  try
    {
      cvec.at(42);
    }
  catch (const std::out_of_range &)
    {
      deallog << "caught correct exception" << std::endl;
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
