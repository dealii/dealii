// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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



// Just output the constraint matrices of the FE_Q element. Test
// introduced when we started to compute them on the fly, rather than
// precomputing them for a number of elements and storing them in a
// table

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>

#include <string>

#include "../tests.h"

#define PRECISION 2



template <int dim>
void
test(const FE_Q<dim> &fe_q)
{
  deallog << fe_q.get_name() << std::endl;

  const FullMatrix<double> &constraints = fe_q.constraints();

  for (unsigned int i = 0; i < constraints.m(); ++i)
    {
      for (unsigned int j = 0; j < constraints.n(); ++j)
        deallog << 100 * constraints(i, j) << ' ';
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

  // no constraints in 1d, but we had the matrices precomputed up to Q4 for 2d
  // and Q2 for 3d
  for (unsigned int degree = 1; degree <= 4; ++degree)
    test<2>(FE_Q<2>(QIterated<1>(QTrapezoid<1>(), degree)));

  // test the standard version (non-equidistant) as well
  test<2>(FE_Q<2>(4));

  for (unsigned int degree = 1; degree <= 2; ++degree)
    test<3>(FE_Q<3>(degree));

  test<3>(FE_Q<3>(4));

  return 0;
}
