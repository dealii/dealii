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


// Test inplace_vector's type traits


#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include "inplace_vector_common.h"

template <typename T, std::size_t N>
void
test_type_traits()
{
  deallog << "is_trivially_copy_constructible_v<T> = "
          << std::is_trivially_copy_constructible_v<T> << std::endl;
  deallog << "is_trivially_copy_constructible_v<inplace_vector<T, " << N
          << ">> = "
          << std::is_trivially_copy_constructible_v<
               std_cxx26::inplace_vector<T, N>> << std::endl;

  deallog << "is_trivially_move_constructible_v<T> = "
          << std::is_trivially_move_constructible_v<T> << std::endl;
  deallog << "is_trivially_move_constructible_v<inplace_vector<T, " << N
          << ">> = "
          << std::is_trivially_move_constructible_v<
               std_cxx26::inplace_vector<T, N>> << std::endl;

  deallog << "is_trivially_copy_assignable_v<T> = "
          << std::is_trivially_copy_assignable_v<T> << std::endl;
  deallog << "is_trivially_copy_assignable_v<inplace_vector<T, " << N << ">> = "
          << std::is_trivially_copy_assignable_v<
               std_cxx26::inplace_vector<T, N>> << std::endl;

  deallog << "is_trivially_move_assignable_v<T> = "
          << std::is_trivially_move_assignable_v<T> << std::endl;
  deallog << "is_trivially_move_assignable_v<inplace_vector<T, " << N << ">> = "
          << std::is_trivially_move_assignable_v<
               std_cxx26::inplace_vector<T, N>> << std::endl;

  deallog << "is_trivially_destructible_v<T> = "
          << std::is_trivially_destructible_v<T> << std::endl;
  deallog << "is_trivially_destructible_v<inplace_vector<T, " << N << ">> = "
          << std::is_trivially_destructible_v<
               std_cxx26::inplace_vector<T, N>> << std::endl;

  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog << std::endl;
  deallog.push("A");
  test_type_traits<A, 1>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("int");
  test_type_traits<int, 1>();
  deallog.pop();

  deallog << std::endl;
  deallog.push("vector<int>");
  test_type_traits<std::vector<int>, 1>();
  deallog.pop();

  deallog.pop();
}
