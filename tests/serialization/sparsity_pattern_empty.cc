// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2022 by the deal.II authors
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


// Check serialization for SparsityPattern. This used to fail for empty
// sparsity patterns with a segfault.

#include <deal.II/lac/sparsity_pattern.h>

#include "../testmatrix.h"
#include "serialization.h"


void
test()
{
  // Try a sparsity pattern that has nonzero size, but is otherwise empty
  const unsigned int N1 = 5;
  SparsityPattern    sp1((N1 - 1) * (N1 - 1), (N1 - 1) * (N1 - 1), 5);
  sp1.compress();

  // Also try one that is completely empty. This is the case that was
  // shown in the original report.
  SparsityPattern sp2;

  {
    SparsityPattern sp3;
    verify(sp1, sp3);
  }

  {
    SparsityPattern sp3;
    verify(sp2, sp3);
  }
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  // Now also try with the usual harness
  test();

  deallog << "OK" << std::endl;
}
