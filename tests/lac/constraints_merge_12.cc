// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Merge two AffineConstraints<double> objects initialized with IndexSets where
// entries have been added in a way that local_lines extends beyond the end of
// the IndexSet

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"



void
merge_check()
{
  deallog << "Checking AffineConstraints<double>::merge with localized lines"
          << std::endl;

  // set local lines to a very large range that surely triggers an error if
  // the implementation is wrong
  IndexSet local_lines(100000000);
  local_lines.add_range(99999800, 100000000);
  local_lines.compress();

  // works correctly
  AffineConstraints<double> c1(local_lines, local_lines),
    c2(local_lines, local_lines);
  for (types::global_dof_index i = 99999800; i < local_lines.size(); ++i)
    if (i % 2 == 1)
      c1.constrain_dof_to_zero(i);

  c2.constrain_dof_to_zero(99999800);
  c2.constrain_dof_to_zero(99999801);
  c2.constrain_dof_to_zero(99999802);

  // now merge the two and print the results
  c2.merge(c1, AffineConstraints<double>::right_object_wins);
  c2.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  merge_check();
}
