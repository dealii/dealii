// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// generate the constraints for a case where there are nodes that have
// a constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then print this set
// of constraints.
//
// we used to get this case wrong (we simply forgot to output this
// node).


#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"



void
test()
{
  AffineConstraints<double> cm;

  // a "regular" constraint
  cm.add_line(1);
  cm.add_entry(1, 2, 42.);

  // a "singular" constraint
  cm.add_line(4);

  cm.print(deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
