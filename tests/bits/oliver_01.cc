// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Oliver found an example, where sparse_matrix_iterator->value=0 didn't work,
// because the iterator->value expects a double on the right hand side, not an
// integer. If the right hand side is zero, it can also be converted to a
// pointer, which leads to an ambiguity. Fix this by having an additional
// operator= in the iterator/reference class

#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"


int
main()
{
  initlog();

  // this test only needs to compile, not run
  if (false)
    {
      SparseMatrix<double>::iterator *i;
      (*i)->value() = (int)0;
    }

  deallog << "OK" << std::endl;

  return 0;
}
