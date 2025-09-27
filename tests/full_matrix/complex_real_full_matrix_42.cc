// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check FullMatrix::invert. like the full_matrix_* tests, but use
// complex-valued matrices and vectors, even though we only store real values
// in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  for (unsigned int n = 1; n <= 5; ++n)
    {
      const std::complex<number> array[] = {50, 2, 3, 4,  5, 6, 50, 8, 9,
                                            0,  1, 2, 50, 4, 5, 6,  7, 8,
                                            50, 0, 1, 2,  3, 4, 50};

      FullMatrix<std::complex<number>> m(n, n, array), p(n, n);
      p.invert(m);
      print_matrix(p);
    }
}
