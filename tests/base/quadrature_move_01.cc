// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
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


template <class Quad, typename... Args>
std::string
check_q_move(Args &&...args)
{
  Quad               quad1(args...);
  const unsigned int size1 = quad1.size();

  std::vector<double> weights1 = quad1.get_weights();

  Quad               quad2(std::move(quad1));
  const unsigned int size2 = quad2.size();

  std::vector<double> weights2 = quad2.get_weights();

  if (size1 != size2)
    return "NOPE";

  for (unsigned short i = 0; i < size1; ++i)
    if (std::fabs(weights1[i] - weights2[i]) > 1.0e-16)
      return "NOPE";

  return "OK";
}


template <template <int dim> class Quad, typename... Args>
void
check_quadrature_move(Args &&...args)
{
  deallog << check_q_move<Quad<1>>(std::forward<Args>(args)...) << 1 << ' '
          << check_q_move<Quad<2>>(std::forward<Args>(args)...) << 2 << ' '
          << check_q_move<Quad<3>>(std::forward<Args>(args)...) << 3
          << std::endl;
}


int
main()
{
  initlog();

  check_quadrature_move<QMidpoint>();
  check_quadrature_move<QTrapezoid>();
  check_quadrature_move<QSimpson>();
  check_quadrature_move<QMilne>();
  check_quadrature_move<QWeddle>();

  for (unsigned int p = 2; p < 5; ++p)
    {
      check_quadrature_move<QGauss>(p);
      check_quadrature_move<QGaussLobatto>(p);
    }

  const auto ep = QGaussRadauChebyshev<1>::right;
  for (unsigned int p = 2; p < 5; ++p)
    {
      deallog << "Gauss Log R: " << check_q_move<QGaussLogR<1>>(p) << std::endl;
      deallog << "Gauss Radau Chebyshev: "
              << check_q_move<QGaussRadauChebyshev<1>>(p, ep) << std::endl;
    }

  return 0;
}
