// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Output is the constraint matrices of the RT_Bubbles element

#include <deal.II/fe/fe_rt_bubbles.h>

#include <string>

#include "../tests.h"

#define PRECISION 8


template <int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_RT_Bubbles<" << dim << "> (" << degree << ')' << std::endl;

  FE_RT_Bubbles<dim>        fe_rt_bubbles(degree);
  const FullMatrix<double> &constraints = fe_rt_bubbles.constraints();

  for (unsigned int i = 0; i < constraints.m(); ++i)
    {
      for (unsigned int j = 0; j < constraints.n(); ++j)
        deallog << constraints(i, j) << ' ';
      deallog << std::endl;
    }

  deallog << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  for (unsigned int degree = 1; degree < 4; ++degree)
    {
      test<2>(degree);
      test<3>(degree);
    }

  return 0;
}
