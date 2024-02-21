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



// Just output the embedding matrices of the FE_DGPMonomial
// element. Test introduced when we started to compute them on the
// fly, rather than precomputing them for a number of elements and
// storing them in a table

#include <deal.II/fe/fe_dgp_monomial.h>

#include <string>

#include "../tests.h"

#define PRECISION 4



template <int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_DGPMonomial<" << dim << "> (" << degree << ')' << std::endl;

  FE_DGPMonomial<dim> fe_q(degree);

  for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell; ++c)
    {
      const FullMatrix<double> &m = fe_q.get_prolongation_matrix(c);

      for (unsigned int i = 0; i < m.m(); ++i)
        {
          for (unsigned int j = 0; j < m.n(); ++j)
            deallog << m(i, j) << ' ';
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

  // we had the matrices precomputed
  // up to Q4 for 1d, Q3 for 2d and
  // Q2 for 3d
  for (unsigned int degree = 1; degree <= 4; ++degree)
    test<1>(degree);

  for (unsigned int degree = 1; degree <= 4; ++degree)
    test<2>(degree);

  for (unsigned int degree = 1; degree <= 2; ++degree)
    test<3>(degree);

  return 0;
}
