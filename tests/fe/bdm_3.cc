// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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



// Just output the constraint matrices of the BDM element

#include <deal.II/fe/fe_bdm.h>

#include <string>

#include "../tests.h"

#define PRECISION 8



template <int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_BDM<" << dim << "> (" << degree << ")" << std::endl;

  FE_BDM<dim>               fe_rt(degree);
  const FullMatrix<double> &constraints = fe_rt.constraints();

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
