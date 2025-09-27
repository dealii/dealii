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


// like_02, but don't just check TableBase::fill but instead use the
// constructors that already initialize


#include <deal.II/base/table.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(0);

  // rank=1
  {
    deallog << "rank=1" << std::endl;

    const double entries[] = {1, 2, 3};

    {
      Table<1, double> t(3, entries, true);

      for (unsigned int i = 0; i < t.size()[0]; ++i)
        deallog << t[i] << ' ';
      deallog << std::endl;
    }

    // passing false as second argument shouldn't
    // make a difference for rank-1 tables
    {
      Table<1, double> t(3, entries, false);
      for (unsigned int i = 0; i < t.size()[0]; ++i)
        deallog << t[i] << ' ';
      deallog << std::endl;
    }
  }


  // rank=2
  {
    deallog << "rank=2" << std::endl;

    const double entries[] = {1, 2, 3, 4, 5, 6};

    // create a 2x3 table from this
    {
      Table<2, double> t(2, 3, entries, true);

      for (unsigned int i = 0; i < t.size()[0]; ++i)
        {
          for (unsigned int j = 0; j < t.size()[1]; ++j)
            deallog << t[i][j] << ' ';
          deallog << std::endl;
        }
    }

    // same data, same table, but filled in transpose ordering
    {
      Table<2, double> t(2, 3, entries, false);

      for (unsigned int i = 0; i < t.size()[0]; ++i)
        {
          for (unsigned int j = 0; j < t.size()[1]; ++j)
            deallog << t[i][j] << ' ';
          deallog << std::endl;
        }
    }
  }


  // rank=3
  {
    deallog << "rank=3" << std::endl;

    const double entries[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    // create a 2x3x2 table from this
    {
      Table<3, double> t(2, 3, 2, entries, true);

      for (unsigned int i = 0; i < t.size()[0]; ++i)
        {
          for (unsigned int j = 0; j < t.size()[1]; ++j)
            {
              deallog << '(';
              for (unsigned int k = 0; k < t.size()[2]; ++k)
                deallog << t[i][j][k] << ' ';
              deallog << ')';
            }
          deallog << std::endl;
        }
    }

    // same data, same table, but filled in transpose ordering
    {
      Table<3, double> t(2, 3, 2, entries, false);

      for (unsigned int i = 0; i < t.size()[0]; ++i)
        {
          for (unsigned int j = 0; j < t.size()[1]; ++j)
            {
              deallog << '(';
              for (unsigned int k = 0; k < t.size()[2]; ++k)
                deallog << t[i][j][k] << ' ';
              deallog << ')';
            }
          deallog << std::endl;
        }
    }
  }
}
