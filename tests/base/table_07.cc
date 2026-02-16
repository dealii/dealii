// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


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
