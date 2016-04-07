// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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


#include "../tests.h"
#include <fstream>
#include <cmath>

#include <deal.II/base/quadrature_lib.h>


template <class Quad, typename... Args>
std::string check_q_move(Args &&...args)
{
  Quad quad1(args...);
  const unsigned int size1 = quad1.size();

  std::vector<double> weights1 = quad1.get_weights();

  Quad quad2(std::move(quad1));
  const unsigned int size2 = quad2.size();

  std::vector<double> weights2 = quad2.get_weights();

  if (size1 != size2) return "NOPE";

  for (unsigned short i = 0; i < size1; ++i)
    if (std::fabs(weights1[i] - weights2[i]) > 1.0e-16)
      return "NOPE";

  return "OK";
}


template <template <int dim> class Quad, typename... Args>
void check_quadrature_move(Args &&...args)
{
  deallog << check_q_move<Quad<1>>(std::forward<Args>(args)...) << 1 << " "
          << check_q_move<Quad<2>>(std::forward<Args>(args)...) << 2 << " "
          << check_q_move<Quad<3>>(std::forward<Args>(args)...) << 3
          << std::endl;
}


int main()
{
  initlog();

  check_quadrature_move<QMidpoint>();
  check_quadrature_move<QTrapez>();
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
