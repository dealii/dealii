// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test AffineConstraints::copy_from() for different types (float, double).

#include <deal.II/lac/affine_constraints.h>

#include <vector>

#include "../tests.h"

template <typename T1, typename T2>
void
test()
{
  AffineConstraints<T1> constraints;

  constraints.add_line(1);
  constraints.add_entry(1, 2, 1.);
  constraints.add_entry(1, 3, 1.);

  constraints.add_line(3);
  constraints.add_entry(3, 4, 1.);
  constraints.add_entry(3, 5, 1.);

  constraints.add_line(5);
  constraints.add_entry(5, 0, 1.);

  constraints.close();

  AffineConstraints<T2> constraints_float;

  constraints_float.copy_from(constraints);

  constraints.print(deallog.get_file_stream());
  deallog << std::endl;
  constraints_float.print(deallog.get_file_stream());
  deallog << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<double, double>();
  test<float, double>();
  test<double, float>();
  test<float, float>();
}
