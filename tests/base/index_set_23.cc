// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
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


// we managed to get a function as simple as IndexSet::operator== wrong -- who
// knew?

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100), is2(100);

  is1.add_range(0, 10);
  is2.add_range(0, 20);

  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is1 != is2) == true, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
