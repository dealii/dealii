// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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

  deallog << "OK" << std::endl;
}
