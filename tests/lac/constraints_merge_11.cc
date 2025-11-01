// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Try to merge two empty AffineConstraints<double> objects initialized with
// IndexSets that don't match.

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"



void
merge_check()
{
  deallog << "Checking AffineConstraints<double>::merge with localized lines"
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
  AffineConstraints<double> c1(local_lines1, local_lines1),
    c2(local_lines2, local_lines2);

  // now merge the two and print the
  // results
  c1.merge(c2, AffineConstraints<double>::no_conflicts_allowed, true);
  c1.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  merge_check();
}
