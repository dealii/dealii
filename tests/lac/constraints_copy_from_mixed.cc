// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
