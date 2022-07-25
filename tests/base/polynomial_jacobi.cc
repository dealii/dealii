// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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

// check jacobi_polynomial_value and jacobi_polynomial_roots

#include <deal.II/base/polynomial.h>

#include "../tests.h"

using namespace Polynomials;


int
main()
{
  initlog();
  deallog.precision(10);

  for (int alpha = 0; alpha < 3; ++alpha)
    for (int beta = 0; beta < 3; ++beta)
      for (unsigned int degree = 0; degree < 40; ++degree)
        {
          deallog << "Jacobi_" << degree << "^(" << alpha << ',' << beta
                  << ") at 0.3: "
                  << jacobi_polynomial_value(degree, alpha, beta, 0.3)
                  << std::endl;

          deallog << "Roots: ";
          std::vector<double> roots =
            jacobi_polynomial_roots<double>(degree, alpha, beta);

          // assert that roots are increasing
          for (unsigned int i = 1; i < roots.size(); ++i)
            Assert(roots[i] > roots[i - 1], ExcInternalError());
          Assert(roots.size() == 0 || roots.front() > -1., ExcInternalError());
          Assert(roots.size() == 0 || roots.back() < 1., ExcInternalError());

          for (unsigned int i = 0; i < roots.size(); ++i)
            deallog << roots[i] << ' ';
          deallog << std::endl << std::endl;
        }
}
