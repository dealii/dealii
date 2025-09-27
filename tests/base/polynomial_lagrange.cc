// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test Lagrange interpolation

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
              if (std::fabs(y - 1.) > 2.e-10)
                deallog << "Error1  lg y=" << std::log10(std::fabs(y - 1.))
                        << std::endl;
            }
          else
            {
              if (std::fabs(y) > 3.e-10)
                deallog << "Error0  lg y=" << std::log10(std::fabs(y))
                        << std::endl;
            }
        }
      deallog << std::endl;
    }
}


void
check_poly(const Quadrature<1> &q)
{
  deallog << "Points: " << q.size() << std::endl;
  std::vector<Polynomial<double>> p =
    generate_complete_Lagrange_basis(q.get_points());
  check_interpolation(p, q.get_points());
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
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  QTrapezoid<1> trapez;
  QSimpson<1>   simpson;
  QIterated<1>  equi7(trapez, 6);
  QIterated<1>  equi10(trapez, 9);

  QGauss<1> g2(2);
  QGauss<1> g3(3);
  QGauss<1> g7(7);
  QGauss<1> g10(10);

  QGaussLobatto<1> gl2(2);
  QGaussLobatto<1> gl3(3);
  QGaussLobatto<1> gl4(4);
  QGaussLobatto<1> gl5(5);
  QGaussLobatto<1> gl6(6);
  QGaussLobatto<1> gl7(7);
  QGaussLobatto<1> gl10(10);

  deallog.push("LagrangeEquidistant");
  check_lge(6);
  check_lge(9);
  deallog.pop();
  deallog.push("Equidistant");
  check_poly(trapez);
  check_poly(simpson);
  check_poly(equi7);
  check_poly(equi10);
  deallog.pop();
  deallog.push("Gauss");
  check_poly(g2);
  check_poly(g3);
  check_poly(g7);
  check_poly(g10);
  deallog.pop();
  deallog.push("GaussLobatto");
  check_poly(gl2);
  check_poly(gl3);
  check_poly(gl4);
  check_poly(gl5);
  check_poly(gl6);
  check_poly(gl7);
  check_poly(gl10);
  deallog.pop();
}
