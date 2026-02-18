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



// check FullMatrix::matrix_scalar_product


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<number> m;
  make_square_matrix(m);
  Vector<number> v, w;
  make_range_vector(v);
  make_range_vector(w);
  for (unsigned int i = 0; i < w.size(); ++i)
    w(i) = w(i) + 1.;

  deallog << m.matrix_scalar_product(v, w) << std::endl;
}
