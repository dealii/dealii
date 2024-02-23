// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
