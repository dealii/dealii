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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Check that pop_back() for `IndexSet` works properly

#include <deal.II/base/index_set.h>

#include "../tests.h"

int
main()
{
  initlog();

  IndexSet is1(10);
  is1.add_range(0, 2);
  is1.add_range(5, 8);

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_back() << std::endl;

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_back() << std::endl;

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_back() << std::endl;


  is1.add_index(9);

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_back() << std::endl;

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_back() << std::endl;
}
