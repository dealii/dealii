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

  // copy the object above
  AffineConstraints<double> constraints_2;
  constraints_2.copy_from(constraints);

  // let both objects describe themselves in string form
  std::ostringstream s1, s2;
  constraints.print(s1);
  constraints_2.print(s2);

  // make sure they're the same
  Assert(s1.str() == s2.str(), ExcInternalError());
}

int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
