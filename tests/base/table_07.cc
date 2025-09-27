// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that Table<N,T> works for types that can not be copy constructed

#include <deal.II/base/table.h>

#include "../tests.h"


class T
{
public:
  T()
  {
    deallog << "Default construct." << std::endl;
  }
  T(const T &) = delete;
  T(T &&)
  {
    deallog << "Move construct." << std::endl;
  }

  T &
  operator=(const T &) = delete;
  T &
  operator=(T &&)
  {
    deallog << "Move assign." << std::endl;
    return *this;
  }
};


int
main()
{
  initlog();
  dealii::Table<2, T> table(2, 2);

  dealii::Table<2, T> table2;
  table2 = std::move(table); // should not create new objects

  table.clear();
  Assert(table.empty() == true, ExcInternalError());

  deallog << "OK" << std::endl;
}
