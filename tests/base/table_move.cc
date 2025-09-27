// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/table.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int     entries[] = {1, 2, 3, 4, 5, 6};
  Table<2, unsigned int> t(2, 3);
  t.fill(entries, true);

  deallog << "Sizes: " << t.size(0) << ", " << t.size(1) << std::endl;
  deallog << "Contents:" << std::endl;
  for (unsigned int i = 0; i < t.size(0); ++i)
    {
      for (unsigned int j = 0; j < t.size(1); ++j)
        deallog << t[i][j] << ' ';
      deallog << std::endl;
    }

  Table<2, unsigned int> s = std::move(t);

  deallog << "Sizes of moved-to table: " << s.size(0) << ", " << s.size(1)
          << std::endl;
  deallog << "Sizes of moved-from table: " << t.size(0) << ", " << t.size(1)
          << std::endl;

  deallog << "Contents of moved-to table:" << std::endl;
  for (unsigned int i = 0; i < s.size(0); ++i)
    {
      for (unsigned int j = 0; j < s.size(1); ++j)
        deallog << s[i][j] << ' ';
      deallog << std::endl;
    }

  const TableIndices<2> new_size(4, 2);
  s.reinit(new_size);
  const unsigned int new_entries[] = {1, 2, 3, 4, 5, 6, 7, 8};
  s.fill(new_entries, true);

  deallog << "Sizes of new table: " << s.size(0) << ", " << s.size(1)
          << std::endl;
  deallog << "Contents of new table:" << std::endl;
  for (unsigned int i = 0; i < s.size(0); ++i)
    {
      for (unsigned int j = 0; j < s.size(1); ++j)
        deallog << s[i][j] << ' ';
      deallog << std::endl;
    }

  t = std::move(s);

  deallog << "Sizes of moved-to table: " << t.size(0) << ", " << t.size(1)
          << std::endl;
  deallog << "Sizes of moved-from table: " << s.size(0) << ", " << s.size(1)
          << std::endl;

  deallog << "Contents of moved-to table:" << std::endl;
  for (unsigned int i = 0; i < t.size(0); ++i)
    {
      for (unsigned int j = 0; j < t.size(1); ++j)
        deallog << t[i][j] << ' ';
      deallog << std::endl;
    }

  return 0;
}
