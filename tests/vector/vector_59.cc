// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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



// check Vector initialization with std::initializer_list objects

#include <deal.II/lac/vector.h>

#include "../tests.h"



int
main()
{
  initlog();

  Vector<double> vd({1.0, 2.0, 3.0});
  for (const auto &x : vd)
    deallog << x << ' ';
  deallog << std::endl;


  Vector<double> vi({1, 2, 3});
  for (const auto &x : vi)
    deallog << x << ' ';
  deallog << std::endl;
}
