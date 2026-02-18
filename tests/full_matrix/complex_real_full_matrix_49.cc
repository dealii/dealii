// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check FullMatrix::precondition_Jacobi. like the full_matrix_* tests, but use
// complex-valued matrices and vectors, even though we only store real values
// in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<std::complex<number>> m;
  make_square_matrix(m);
  Vector<std::complex<number>> v, w;
  make_range_vector(v);
  make_range_vector(w);

  m.precondition_Jacobi(v, w, 3.141);
  print_vector(v);
  print_vector(w);
}
