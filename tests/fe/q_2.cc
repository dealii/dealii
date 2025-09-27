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



// Just output the embedding matrices of the FE_Q element. Test
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

  for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell; ++c)
    {
      const FullMatrix<double> &m = fe_q.get_prolongation_matrix(c);

      for (unsigned int i = 0; i < m.m(); ++i)
        {
          for (unsigned int j = 0; j < m.n(); ++j)
            deallog << 100 * m(i, j) << ' ';
          deallog << std::endl;
        }

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

  // we had the matrices precomputed up to Q4 for 1d, Q3 for 2d and Q2 for 3d
  for (unsigned int degree = 1; degree <= 4; ++degree)
    test<1>(FE_Q<1>(QIterated<1>(QTrapezoid<1>(), degree)));

  // test the standard version (non-equidistant) as well
  test<1>(FE_Q<1>(4));

  for (unsigned int degree = 1; degree <= 3; ++degree)
    test<2>(FE_Q<2>(QIterated<1>(QTrapezoid<1>(), degree)));

  test<2>(FE_Q<2>(4));

  for (unsigned int degree = 1; degree <= 2; ++degree)
    test<3>(FE_Q<3>(degree));

  test<3>(FE_Q<3>(4));

  return 0;
}
