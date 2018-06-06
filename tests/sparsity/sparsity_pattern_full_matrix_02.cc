
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// Check SparsityPattern::copy_from(FullMatrix). The peer review of
// the original patch uncovered a case where we would segfault with
// empty matrices. Cover this in a testcase as well.

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


void
test(const unsigned int M, const unsigned int N)
{
  SparsityPattern sp(M, N);

  // fill a full matrix completely
  FullMatrix<double> mat(M, N);
  for (unsigned int m = 0; m < M; ++m)
    for (unsigned int n = 0; n < N; ++n)
      mat(m, n) = 1;

  // then copy the nonzero entries of 'mat' (i.e., all entries) to the
  // sparsity pattern
  sp.copy_from(mat);

  deallog << sp.n_nonzero_elements() << std::endl;
}



int
main()
{
  initlog();

  test(0, 0);
}
