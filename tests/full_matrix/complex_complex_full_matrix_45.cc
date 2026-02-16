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



// check FullMatrix::vmult. like the full_matrix_* tests, but use
// complex-valued matrices and vectors; this time we actually store complex
// values in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<std::complex<number>> m;
  make_complex_matrix(m);
  Vector<std::complex<number>> v, w;
  make_complex_range_vector(v);
  make_complex_domain_vector(w);

  m.vmult(v, w, true);
  print_vector(v);
  print_vector(w);

  m.vmult(v, w, false);
  print_vector(v);
  print_vector(w);
}
