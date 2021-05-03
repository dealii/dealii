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



// generate the constraints for a case where there are nodes that have
// a constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then check that if
// we condense a matrix with this then the right thing happens


#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"



void
test()
{
  // constrain each dof to zero. this
  // should yield a diagonal matrix
  AffineConstraints<double> cm;

  for (unsigned int i = 0; i < 5; ++i)
    cm.add_line(i);
  cm.close();

  // completely fill a 5x5 matrix
  // with some values
  SparsityPattern sp(5, 5, 5);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      sp.add(i, j);
  sp.compress();

  SparseMatrix<double> m(sp);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      m.set(i, j, i + j + 2);

  // now condense it
  cm.condense(m);

  // and print it
  for (unsigned int i = 0; i < 5; ++i)
    {
      for (unsigned int j = 0; j < 5; ++j)
        deallog << m(i, j) << ' ';
      deallog << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
