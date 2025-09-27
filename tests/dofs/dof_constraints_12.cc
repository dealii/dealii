// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that the AffineConstraints::shift() method works as expected.

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


int
main()
{
  initlog();

  unsigned int n1 = 10, n2 = 5;

  auto local_lines1 = complete_index_set(n1);
  auto local_lines2 = complete_index_set(n2);
  auto local_lines  = complete_index_set(n1 + n2);

  // Fake two independent constraints on two independent dofs.
  AffineConstraints<double> constraints1; //(local_lines1);
  AffineConstraints<double> constraints2; //(local_lines2);
  AffineConstraints<double> constraints;  //(local_lines);

  constraints1.add_line(1);
  constraints1.add_entry(1, 2, 1.0);

  constraints2.add_line(2);
  constraints2.add_entry(2, 3, 2.0);

  constraints1.close();
  constraints2.close();

  deallog << "constraints1: " << std::endl;
  constraints1.print(deallog.get_file_stream());
  deallog << "constraints2: " << std::endl;
  constraints2.print(deallog.get_file_stream());

  // Now I want to build the union of the two constraints.
  // According to the documentation, I can shift the second, then merge.
  // Since in applications you usually need also the second constraints,
  // do this on a copy.
  {
    AffineConstraints<double> tmp; //(local_lines2);
    tmp.merge(constraints2);
    tmp.shift(n1);
    deallog << "constraints2 shifted:" << std::endl;
    tmp.print(deallog.get_file_stream());

    constraints.merge(constraints1,
                      AffineConstraints<double>::no_conflicts_allowed,
                      true);
    constraints.merge(tmp,
                      AffineConstraints<double>::no_conflicts_allowed,
                      true);
  }
  constraints.close();
  deallog << "constraints: " << std::endl;
  constraints.print(deallog.get_file_stream());
}
