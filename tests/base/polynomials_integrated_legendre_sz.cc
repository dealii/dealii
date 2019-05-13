// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// This tests the stability of the polynomial evaluation of
// IntegratedLegendreSZ.

#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_integrated_legendre_sz.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


using namespace Polynomials;


void
check_at_one(const std::vector<Polynomial<double>> &p)
{
  // Ignore first two polynomials as the integrated Legendre polynomials are
  // only defined for degree > 1, it is only added to maintain the recursive
  // relation.

  deallog << "Function value of polynomial at right end point: ";
  for (unsigned int i = 2; i < p.size(); ++i)
    {
      deallog << '.';
      const double y = p[i].value(1.);
      if (std::fabs(y) > 1e-13)
        deallog << "Error1  lg y=" << std::log10(std::fabs(y)) << std::endl;
    }
  deallog << std::endl;
}



void
check_at_half(const std::vector<Polynomial<double>> &p)
{
  // Ignore first two polynomials as the integrated Legendre polynomials are
  // only defined for degree > 1, it is only added to maintain the recursive
  // relation.

  deallog << "Function value of polynomial at -0.5 | 0.5:" << std::endl;
  for (unsigned int i = 2; i < p.size(); ++i)
    {
      const double y = p[i].value(0.5);
      const double z = p[i].value(-0.5);
      deallog << y << " | " << z << std::endl;
    }
  deallog << std::endl;
}



void
check_poly(const unsigned int n)
{
  deallog << "Degree: " << n + 1 << std::endl;
  std::vector<Polynomial<double>> p =
    IntegratedLegendreSZ::generate_complete_basis(n);
  check_at_one(p);
  check_at_half(p);
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check_poly(25);
}
