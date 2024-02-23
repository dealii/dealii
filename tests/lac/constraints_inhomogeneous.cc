// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test: set inhomogeneous constraints and indirectly apply those to other
// constrained nodes.


#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"



void
test()
{
  AffineConstraints<double> cm;

  // an inhomogeneous constraint
  cm.add_line(4);
  cm.set_inhomogeneity(4, 3.14159);

  // a homogeneous constraint that is
  // constrained to the inhomogeneous one
  cm.add_line(1);
  cm.add_entry(1, 2, 42.);
  cm.add_entry(1, 4, 1.);

  // and a standard homogeneous constraint
  cm.add_line(17);
  cm.add_entry(17, 6, 2.);
  cm.add_entry(17, 15, 3.);

  // a "singular" constraint
  cm.add_line(3);

  // now close the constraint matrix
  cm.close();

  cm.print(deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
