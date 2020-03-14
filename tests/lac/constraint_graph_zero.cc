// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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



// generate constraints for a case where there are nodes that have a
// constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then write this as
// a DOT graph.
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

  cm.write_dot(deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
