// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check FullMatrix::trace


#include <deal.II/lac/full_matrix.h>

#include "../tests.h"



int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);

  const unsigned int N = 20;
  FullMatrix<double> m(N, N);

  double tr = 0;
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        m(i, j) = i + j;
        if (i == j)
          tr += i + j;
      }

  deallog << "Trace=" << m.trace() << std::endl;
  Assert(m.trace() == tr, ExcInternalError());
}
