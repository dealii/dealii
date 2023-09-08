// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


// Ensure that we can compare IndexSet objects against the empty
// object. (Comparisons between IndexSet objects of different sizes
// are otherwise not allowed.)

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);
  IndexSet is2;

  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is2 != is1) == true, ExcInternalError());
  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is2 != is1) == true, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
