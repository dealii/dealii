// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

// Check that the move constructor for `IndexSet` works properly

#include "../tests.h"
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/index_set.h>

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  IndexSet is1(100);
  is1.add_range(0, 10);
  is1.add_range(30, 40);

  deallog << is1.size() << ", " << is1.n_elements() << std::endl;

  // Test that move construction works correctly and that the moved object is
  // restored to the default state
  IndexSet is2 = std::move(is1);

  deallog << is2.size() << ", " << is2.n_elements() << std::endl;
  deallog << is1.size() << ", " << is1.n_elements() << std::endl;

  // Test that re-initializing the moved variable works
  is1.set_size(200);
  is1.add_range(90, 110);
  is1.add_range(130, 140);
  is1.add_range(145, 150);

  deallog << is1.size() << ", " << is1.n_elements() << std::endl;

  // Test that move assignment works correctly and that the moved object is
  // restored to the default state
  is2 = std::move(is1);
  deallog << is2.size() << ", " << is2.n_elements() << std::endl;
  deallog << is1.size() << ", " << is1.n_elements() << std::endl;
}
