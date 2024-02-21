// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check usage of TransposeTable iterators with the table itself

#include <deal.II/base/table.h>

#include <numeric>

#include "../tests.h"

int
main()
{
  initlog();

  deallog << std::boolalpha;

  // test a non-empty rectangular table
  TransposeTable<double> table(3, 4);
  std::iota(table.begin(), table.end(), 1.0);
  for (const auto &entry : table)
    {
      deallog << entry.row() << ", " << entry.column() << ", " << entry.value()
              << std::endl;
    }

  deallog << "backwards order:" << std::endl;
  auto it = table.end() - 1;
  for (; it >= table.begin(); --it)
    {
      deallog << it->row() << ", " << it->column() << ", " << it->value()
              << std::endl;
    }
  deallog << "iterator is one before the beginning: "
          << (it == table.begin() - 1) << std::endl;

  deallog << "every other entry:" << std::endl;
  it = table.begin();
  for (; it < table.end(); it += 2)
    {
      deallog << it->row() << ", " << it->column() << ", " << it->value()
              << std::endl;
    }

  // print every other entry
  it = table.end() - 1;
  deallog << "every other entry:" << std::endl;
  for (; it >= table.begin(); it -= 2)
    {
      deallog << it->row() << ", " << it->column() << ", " << it->value()
              << std::endl;
    }

  // test some type equalities
  static_assert(
    std::is_same_v<decltype(table.begin()->value()), double &>,
    "The iterator value for a non-const table should not be const.");
  static_assert(
    std::is_same_v<decltype(table.end()->value()), double &>,
    "The iterator value for a non-const table should not be const.");

  const TransposeTable<double> &ref = table;
  static_assert(std::is_same_v<decltype(ref.begin()->value()), const double &>,
                "The iterator value for a constant table should be const.");
  static_assert(std::is_same_v<decltype(ref.end()->value()), const double &>,
                "The iterator value for a constant table should be const.");

  deallog << "OK" << std::endl;
}
