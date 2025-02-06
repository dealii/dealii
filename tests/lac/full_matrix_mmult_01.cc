// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that an assertion is thrown if output argument of FullMatrix::mmult() is
// also an input argument

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  try
    {
      dealii::FullMatrix<double> A(2, 2);
      dealii::FullMatrix<double> B(2, 2);

      // This should trigger an assertion
      A.mmult(A, B);
    }
  catch (const dealii::ExceptionBase &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
    }

  return 0;
}
