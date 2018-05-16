
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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Check SparsityPattern::copy_from(FullMatrix). This function took
// extraordinary amounts of CPU time: for a 1 x 400,000 matrix, it
// took around 50 seconds; after the patch that fixed the quadratic
// behavior, it takes less than one second
//
// To make sure we really catch the old bevahior through the test
// suite time out, we also do a matrix of size 1 x 4,000,000, which
// would have taken more than an hour before and now still takes less
// than one second. Finally, also try a square matrix because this
// case is treated specially in the function.

#include "../tests.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>


void
test (const unsigned int M,
      const unsigned int N)
{
  SparsityPattern sp (M,N);

  // fill a full matrix completely
  FullMatrix<double> mat(M,N);
  for (unsigned int m=0; m<M; ++m)
    for (unsigned int n=0; n<N; ++n)
      mat(m,n) = 1;

  // then copy the nonzero entries of 'mat' (i.e., all entries) to the
  // sparsity pattern
  sp.copy_from (mat);

  deallog << sp.n_nonzero_elements() << std::endl;
}



int
main ()
{
  initlog();

  test (1, 400000);
  test (1, 4000000);
  test (10, 10);
}
