// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Just output the constraint matrices of the FE_Q element. Test
// introduced when we started to compute them on the fly, rather than
// precomputing them for a number of elements and storing them in a
// table

#include <deal.II/fe/fe_q.h>

#include <string>

#include "../tests.h"

#define PRECISION 2



template <int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_Q<" << dim << "> (" << degree << ')' << std::endl;

  FE_Q<dim> fe_q(degree);
  FE_Q<dim> x(fe_q);
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  // no constraints in 1d, but we had
  // the matrices precomputed up to
  // Q4 for 2d and Q2 for 3d
  for (unsigned int degree = 1; degree <= 4; ++degree)
    test<2>(degree);

  for (unsigned int degree = 1; degree <= 2; ++degree)
    test<3>(degree);

  return 0;
}
