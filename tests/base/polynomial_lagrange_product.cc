// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests correctness of values and derivatives for polynomials derived from
// Lagrange product form

#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


using namespace Polynomials;


void
check_derivatives(const std::vector<Polynomial<double>> &p,
                  const unsigned int                     n_deriv)
{
  // check whether the values and derivatives
  // are evaluated correctly some randomly
  // generated points. compare with a polynomial
  // that is not in product form (we get the
  // expanded form by adding a dummy polynomial;
  // addition of polynomials destroys the
  // product form in the current implementation)
  deallog << "Representation of derivatives up to order " << n_deriv - 1
          << std::endl;
  std::vector<double> values(n_deriv), values_ref(n_deriv);
  Monomial<double>    zero(0, 0);
  for (unsigned int j = 0; j < p.size(); ++j)
    {
      double x = random_value<double>();
      p[j].value(x, values);
      Polynomial<double> q = p[j];
      q += zero;
      q.value(x, values_ref);
      for (unsigned int i = 0; i < n_deriv; ++i)
        {
          deallog << '.';
          if (std::fabs(values[i] - values_ref[i]) >
              std::max(1e-11, 1e-11 * std::fabs(values[i])))
            deallog << "Error deriv" << i << "  lg y="
                    << std::log10(std::fabs(values[i] - values_ref[i]))
                    << ", is: " << values[i] << ", should be: " << values_ref[i]
                    << std::endl;
        }
    }
  deallog << std::endl;
}



void
check_poly(const Quadrature<1> &q)
{
  deallog << "Points: " << q.size() << std::endl;
  std::vector<Polynomial<double>> p =
    generate_complete_Lagrange_basis(q.get_points());
  for (unsigned int i = 1; i < 6; ++i)
    check_derivatives(p, i);
  deallog << std::endl;
}



void
check_lge(unsigned int n)
{
  deallog << "Points: " << n + 1 << std::endl;
  std::vector<Polynomial<double>> p =
    LagrangeEquidistant::generate_complete_basis(n);
  for (unsigned int i = 1; i < 6; ++i)
    check_derivatives(p, i);
  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog.push("LagrangeEquidistant");
  for (unsigned i = 1; i < 8; i += 2)
    check_lge(i);
  deallog.pop();
  deallog << std::endl;

  // Lagrange elements on GL points have good
  // conditioning, so test to some very high
  // orders
  deallog.push("GaussLobatto");
  for (unsigned i = 1; i < 8; i += 2)
    check_poly(QGaussLobatto<1>(i + 1));
  deallog.pop();
}
