// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Output is the constraint matrices of the RT_Bubbles element

#include "../tests.h"
#include <deal.II/fe/fe_rt_bubbles.h>

#include <string>

#define PRECISION 8

template <int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_RT_Bubbles<" << dim << "> (" << degree << ")" << std::endl;

  FE_RT_Bubbles<dim>        fe_rt_bubbles(degree);
  const FullMatrix<double>& constraints = fe_rt_bubbles.constraints();

  for(unsigned int i = 0; i < constraints.m(); ++i)
    {
      for(unsigned int j = 0; j < constraints.n(); ++j)
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

  for(unsigned int degree = 1; degree < 4; ++degree)
    {
      test<2>(degree);
      test<3>(degree);
    }

  return 0;
}
