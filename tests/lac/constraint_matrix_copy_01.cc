// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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
