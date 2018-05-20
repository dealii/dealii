// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// Similar to polynomial_lagrange_order, but for Legendre interpolation
// This tests the stability of the polynomial evaluation

#include "../tests.h"

#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

using namespace Polynomials;

void
check_at_one(const std::vector<Polynomial<double>>& p)
{
  deallog << "Function value of polynomial at right end point: ";
  for(unsigned int i = 0; i < p.size(); ++i)
    {
      deallog << '.';
      const double y = p[i].value(1.);
      if(std::fabs(y - std::sqrt(2 * i + 1)) > 1e-13 * std::sqrt(2 * i + 1))
        deallog << "Error1  lg y=" << std::log10(std::fabs(y - 1.))
                << std::endl;
    }
  deallog << std::endl;
}

void
check_poly(const unsigned int n)
{
  deallog << "Degree: " << n + 1 << std::endl;
  std::vector<Polynomial<double>> p = Legendre::generate_complete_basis(n);
  check_at_one(p);
  deallog << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  check_poly(10);
  check_poly(50);
}
