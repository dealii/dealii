// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/affine_constraints.h>

#include <vector>

#include "../tests.h"

void
test()
{
  AffineConstraints<double> constraints;

  constraints.add_line(1);
  constraints.add_entry(1, 2, 1.);
  constraints.add_entry(1, 3, 1.);

  constraints.add_line(3);
  constraints.add_entry(3, 4, 1.);
  constraints.add_entry(3, 5, 1.);

  constraints.add_line(5);
  constraints.add_entry(5, 0, 1.);

  constraints.close();

  std::vector<types::global_dof_index> indices(4);
  indices[0] = 1;
  indices[1] = 2;
  indices[2] = 3;
  indices[3] = 5;

  constraints.resolve_indices(indices);

  for (unsigned int i = 0; i < indices.size(); ++i)
    deallog << "Index: " << indices[i] << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
