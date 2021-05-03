// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
    std::is_same<decltype(table.begin()->value()), double &>::value,
    "The iterator value for a non-const table should not be const.");
  static_assert(
    std::is_same<decltype(table.end()->value()), double &>::value,
    "The iterator value for a non-const table should not be const.");

  const TransposeTable<double> &ref = table;
  static_assert(
    std::is_same<decltype(ref.begin()->value()), const double &>::value,
    "The iterator value for a constant table should be const.");
  static_assert(
    std::is_same<decltype(ref.end()->value()), const double &>::value,
    "The iterator value for a constant table should be const.");

  deallog << "OK" << std::endl;
}
