// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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



// check FullMatrix::invert. like the full_matrix_* tests, but use
// complex-valued matrices and vectors; this time we actually store complex
// values in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  for (unsigned int n = 1; n <= 5; ++n)
    {
      const std::complex<number> array[] = {
        std::complex<number>(50.0, 50.0), std::complex<number>(2.0, 6.0),
        std::complex<number>(3.0, 1.0),   std::complex<number>(4.0, 6.0),
        std::complex<number>(5.0, 1.0),   std::complex<number>(6.0, 2.0),
        std::complex<number>(50.0, 50.0), std::complex<number>(8.0, 2.0),
        std::complex<number>(9.0, 7.0),   std::complex<number>(0.0, 2.0),
        std::complex<number>(1.0, 3.0),   std::complex<number>(2.0, 8.0),
        std::complex<number>(50.0, 50.0), std::complex<number>(4.0, 8.0),
        std::complex<number>(5.0, 3.0),   std::complex<number>(6.0, 4.0),
        std::complex<number>(7.0, 9.0),   std::complex<number>(8.0, 4.0),
        std::complex<number>(50.0, 50.0), std::complex<number>(0.0, 4.0),
        std::complex<number>(1.0, 5.0),   std::complex<number>(2.0, 0.0),
        std::complex<number>(3.0, 5.0),   std::complex<number>(4.0, 0.0),
        std::complex<number>(50.0, 50.0)};

      FullMatrix<std::complex<number>> m(n, n, array), p(n, n);
      p.invert(m);
      print_matrix(p);
    }
}
