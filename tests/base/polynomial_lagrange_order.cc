// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Similar to polynomial_lagrange, but test Lagrange interpolation for high
// order with tighter tolerances, in particular the effect of stability of the
// polynomial evaluation at random points

#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


using namespace Polynomials;


void
check_interpolation(const std::vector<Polynomial<double>> &p,
                    const std::vector<Point<1>>           &x)
{
  for (unsigned int i = 0; i < p.size(); ++i)
    {
      deallog << i;
      for (unsigned int k = 0; k < x.size(); ++k)
        {
          deallog << '.';
          const double y = p[i].value(x[k][0]);
          if (i == k)
            {
              if (std::fabs(y - 1.) > 1e-13)
                deallog << "Error1  lg y=" << std::log10(std::fabs(y - 1.))
                        << std::endl;
            }
          else
            {
              if (std::fabs(y) > 1e-13)
                deallog << "Error0  lg y=" << std::log10(std::fabs(y))
                        << std::endl;
            }
        }
      deallog << std::endl;
    }
}



void
check_constant(const std::vector<Polynomial<double>> &p)
{
  // check whether the sum of all polynomials in
  // the basis gives one for a given point
  deallog << "Representation of one at random points";
  for (unsigned int j = 0; j < 12; ++j)
    {
      double x     = random_value<double>();
      double value = 0;
      for (unsigned int i = 0; i < p.size(); ++i)
        {
          value += p[i].value(x);
        }
      deallog << '.';
      if (std::fabs(1. - value) > 1e-13)
        deallog << "Error1  lg y=" << std::log10(std::fabs(1. - value))
                << std::endl;
    }
  deallog << std::endl;
}



void
check_poly(const Quadrature<1> &q)
{
  deallog << "Points: " << q.size() << std::endl;
  std::vector<Polynomial<double>> p =
    generate_complete_Lagrange_basis(q.get_points());
  check_interpolation(p, q.get_points());
  check_constant(p);
  deallog << std::endl;
}



void
check_lge(unsigned int n)
{
  deallog << "Points: " << n + 1 << std::endl;
  std::vector<Polynomial<double>> p =
    LagrangeEquidistant::generate_complete_basis(n);
  std::vector<Point<1>> x(n + 1);
  const double          h = 1. / n;
  for (unsigned int i = 0; i <= n; ++i)
    x[i][0] = h * i;
  check_interpolation(p, x);
  check_constant(p);
  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog.push("LagrangeEquidistant");
  for (unsigned i = 8; i < 18; i += 2)
    check_lge(i);
  deallog.pop();
  deallog << std::endl;

  // Lagrange elements on GL points have good
  // conditioning, so test to some very high
  // orders
  deallog.push("GaussLobatto");
  for (unsigned i = 8; i < 40; i += 3)
    check_poly(QGaussLobatto<1>(i + 1));
  deallog.pop();
}
