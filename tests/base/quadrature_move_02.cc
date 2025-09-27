// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


template <template <int dim> class Quad, int dim, typename... Args>
std::string
check_q_assign_move(Args &&...args)
{
  Quad<dim>                     quad1(args...);
  const unsigned int            size1    = quad1.size();
  const std::vector<double>     weights1 = quad1.get_weights();
  const std::vector<Point<dim>> points1  = quad1.get_points();

  Quadrature<dim> quad2;
  AssertThrow(quad2.size() == 0, ExcInternalError());

  quad2 = std::move(quad1);

  AssertThrow(quad1.size() == 0, ExcInternalError());
  AssertThrow(quad2.size() == size1, ExcInternalError());

  const std::vector<double>     weights2 = quad2.get_weights();
  const std::vector<Point<dim>> points2  = quad2.get_points();
  for (unsigned int i = 0; i < size1; ++i)
    {
      AssertThrow(std::abs(weights1[i] - weights2[i]) < 1.e-16,
                  ExcInternalError());
      AssertThrow((points1[i] - points2[i]).norm() < 1.e-16,
                  ExcInternalError());
    }

  return "OK";
}


template <template <int dim> class Quad, typename... Args>
void
check_quadrature_assign_move(Args &&...args)
{
  deallog << check_q_assign_move<Quad, 1>(std::forward<Args>(args)...) << 1
          << ' ' << check_q_assign_move<Quad, 2>(std::forward<Args>(args)...)
          << 2 << ' '
          << check_q_assign_move<Quad, 3>(std::forward<Args>(args)...) << 3
          << std::endl;
}


int
main()
{
  initlog();

  check_quadrature_assign_move<QMidpoint>();
  check_quadrature_assign_move<QTrapezoid>();
  check_quadrature_assign_move<QSimpson>();
  check_quadrature_assign_move<QMilne>();
  check_quadrature_assign_move<QWeddle>();

  for (unsigned int p = 2; p < 5; ++p)
    {
      check_quadrature_assign_move<QGauss>(p);
      check_quadrature_assign_move<QGaussLobatto>(p);
    }

  const auto ep = QGaussRadauChebyshev<1>::right;
  for (unsigned int p = 2; p < 5; ++p)
    {
      deallog << "Gauss Log R: " << check_q_assign_move<QGaussLogR, 1>(p)
              << std::endl;
      deallog << "Gauss Radau Chebyshev: "
              << check_q_assign_move<QGaussRadauChebyshev, 1>(p, ep)
              << std::endl;
    }

  return 0;
}
