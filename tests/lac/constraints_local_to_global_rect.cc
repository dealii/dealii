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



// this function tests the correctness of the implementation of
// AffineConstraints<double>::distribute_local_to_global for row and column
// indices with different constraints on the rows and columns.

#include <deal.II/base/function.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>

#include <complex>
#include <iostream>

#include "../tests.h"

void
test()
{
  FullMatrix<double> local(2, 2);
  local(0, 0) = 8.;
  local(0, 1) = 2.;
  local(1, 0) = -2.;
  local(1, 1) = 12.;
  FullMatrix<double>        global1(4, 6), global2(4, 6);
  AffineConstraints<double> cm1, cm2;
  cm1.add_line(2);
  cm1.add_entry(2, 1, 0.5);
  cm1.add_entry(2, 3, 0.5);
  cm1.close();
  cm2.add_line(4);
  cm2.add_entry(4, 2, 0.3);
  cm2.add_entry(4, 5, 0.7);
  cm2.close();

  std::vector<types::global_dof_index> indices1(2), indices2(2);
  indices1[0] = 1;
  indices1[1] = 2;
  indices2[0] = 4;
  indices2[1] = 5;
  cm1.distribute_local_to_global(local, indices1, indices2, global1);
  deallog << "Same constraints: " << std::endl;
  global1.print_formatted(deallog.get_file_stream(), 2, true, 0, "0");
  cm1.distribute_local_to_global(local, indices1, cm2, indices2, global2);
  deallog << "Different constraints: " << std::endl;
  global2.print_formatted(deallog.get_file_stream(), 2, true, 0, "0");
}


int
main()
{
  initlog();

  test();
}
