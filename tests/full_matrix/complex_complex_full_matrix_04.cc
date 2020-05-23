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



// check creation and output of a matrix using an array. like the full_matrix_*
// tests, but use complex-valued matrices and vectors; this time we actually
// store complex values in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  const std::complex<number> array[] = {
    std::complex<number>(1., 1. + 1), std::complex<number>(2., 2. + 1),
    std::complex<number>(3., 3. + 1), std::complex<number>(4., 4. + 1),
    std::complex<number>(5., 5. + 1),

    std::complex<number>(6., 6. + 1), std::complex<number>(7., 7. + 1),
    std::complex<number>(8., 8. + 1), std::complex<number>(9., 9. + 1),
    std::complex<number>(0., 0. + 1),

    std::complex<number>(1., 1. + 1), std::complex<number>(2., 2. + 1),
    std::complex<number>(3., 3. + 1), std::complex<number>(4., 4. + 1),
    std::complex<number>(5., 5. + 1),

    std::complex<number>(6., 6. + 1), std::complex<number>(7., 7. + 1),
    std::complex<number>(8., 8. + 1), std::complex<number>(9., 9. + 1),
    std::complex<number>(0., 0. + 1),

    std::complex<number>(1., 1. + 1), std::complex<number>(2., 2. + 1),
    std::complex<number>(3., 3. + 1), std::complex<number>(4., 4. + 1),
    std::complex<number>(5., 5. + 1)};

  FullMatrix<std::complex<number>> m(5, 5, array);

  print_matrix(m);
}
