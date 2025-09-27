// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for SparsityPattern

#include <deal.II/lac/sparsity_pattern.h>

#include "../testmatrix.h"
#include "serialization.h"


void
test()
{
  const unsigned int N1 = 5;
  SparsityPattern    sp1((N1 - 1) * (N1 - 1), (N1 - 1) * (N1 - 1), 5);
  FDMatrix(N1, N1).five_point_structure(sp1);
  sp1.compress();

  const unsigned int N2 = 3;
  SparsityPattern    sp2((N2 - 1) * (N2 - 1), (N2 - 1) * (N2 - 1), 5);
  FDMatrix(N2, N2).five_point_structure(sp2);
  sp2.compress();

  SparsityPattern sp3;

  verify(sp1, sp2);

  verify(sp1, sp3);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
