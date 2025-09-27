// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// previously, the FETools::get_interpolation_matrix function would
// compute its result itself by interpolation. now, the different
// finite elements do that themselves, if they can. make sure the
// result doesn't change

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <string>

#include "../tests.h"

#define PRECISION 5



template <int dim>
void
test(const FE_Q<dim> fe1, const FE_Q<dim> fe2)
{
  deallog << fe1.get_name() << " to " << fe2.get_name() << std::endl;

  FullMatrix<float> m(fe2.dofs_per_cell, fe1.dofs_per_cell);
  FETools::get_interpolation_matrix(fe1, fe2, m);

  for (unsigned int i = 0; i < m.m(); ++i)
    {
      for (unsigned int j = 0; j < m.n(); ++j)
        deallog << m(i, j) << ' ';

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

  for (unsigned int degree1 = 1; degree1 <= 4; ++degree1)
    for (unsigned int degree2 = 1; degree2 <= 4; ++degree2)
      test<1>(FE_Q<1>(QIterated<1>(QTrapezoid<1>(), degree1)),
              FE_Q<1>(QIterated<1>(QTrapezoid<1>(), degree2)));
  for (unsigned int degree1 = 2; degree1 <= 5; ++degree1)
    for (unsigned int degree2 = 2; degree2 <= 5; ++degree2)
      test<1>(FE_Q<1>(QGaussLobatto<1>(degree1)),
              FE_Q<1>(QGaussLobatto<1>(degree2)));

  for (unsigned int degree1 = 1; degree1 <= 3; ++degree1)
    for (unsigned int degree2 = 1; degree2 <= 3; ++degree2)
      test<2>(FE_Q<2>(QIterated<1>(QTrapezoid<1>(), degree1)),
              FE_Q<2>(QIterated<1>(QTrapezoid<1>(), degree2)));
  for (unsigned int degree1 = 2; degree1 <= 4; ++degree1)
    for (unsigned int degree2 = 2; degree2 <= 4; ++degree2)
      test<2>(FE_Q<2>(QGaussLobatto<1>(degree1)),
              FE_Q<2>(QGaussLobatto<1>(degree2)));

  for (unsigned int degree1 = 1; degree1 <= 2; ++degree1)
    for (unsigned int degree2 = 1; degree2 <= 2; ++degree2)
      test<3>(FE_Q<3>(QIterated<1>(QTrapezoid<1>(), degree1)),
              FE_Q<3>(QIterated<1>(QTrapezoid<1>(), degree2)));
  for (unsigned int degree1 = 2; degree1 <= 4; ++degree1)
    for (unsigned int degree2 = 2; degree2 <= 4; ++degree2)
      test<3>(FE_Q<3>(QGaussLobatto<1>(degree1)),
              FE_Q<3>(QGaussLobatto<1>(degree2)));


  return 0;
}
