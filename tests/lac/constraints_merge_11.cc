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

// Try to merge two empty ConstraintMatrix objects initialized with IndexSets
// that don't match.

#include "../tests.h"
#include <deal.II/lac/constraint_matrix.h>

void
merge_check()
{
  deallog << "Checking ConstraintMatrix::merge with localized lines"
          << std::endl;

  // set local lines to a very large range that
  // surely triggers an error if the
  // implementation is wrong
  IndexSet local_lines1(100000000);
  local_lines1.add_range(99999890, 99999900);
  local_lines1.add_range(99999990, 99999992);
  local_lines1.compress();
  local_lines1.print(deallog.get_file_stream());

  IndexSet local_lines2(100000000);
  local_lines2.add_range(99999893, 99999900);
  local_lines2.add_range(99999990, 100000000);
  local_lines2.compress();
  local_lines2.print(deallog.get_file_stream());

  // works correctly
  ConstraintMatrix c1(local_lines1), c2(local_lines2);

  // now merge the two and print the
  // results
  c1.merge(c2, ConstraintMatrix::no_conflicts_allowed, true);
  c1.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  merge_check();
}
