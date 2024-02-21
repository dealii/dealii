// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check TableBase::fill using an istream_iterator


#include <deal.II/base/table.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(0);

  const std::string elements = "1 2 3 4 5 6";
  {
    // create a 2x3 table from this
    Table<2, double>   t(2, 3);
    std::istringstream in1(elements);
    t.fill(std::istream_iterator<double>(in1), true);

    for (unsigned int i = 0; i < t.size()[0]; ++i)
      {
        for (unsigned int j = 0; j < t.size()[1]; ++j)
          deallog << t[i][j] << ' ';
        deallog << std::endl;
      }

    // same data, same table, but filled in transpose ordering
    std::istringstream in2(elements);
    t.fill(std::istream_iterator<double>(in2), false);

    for (unsigned int i = 0; i < t.size()[0]; ++i)
      {
        for (unsigned int j = 0; j < t.size()[1]; ++j)
          deallog << t[i][j] << ' ';
        deallog << std::endl;
      }
  }
}
