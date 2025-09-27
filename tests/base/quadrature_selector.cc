// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// make sure that the QuadratureSelector works for a selection of
// arguments


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>

#include <string>

#include "../tests.h"


template <int dim>
void
check(const std::string     &name,
      const unsigned int     order,
      const Quadrature<dim> &q)
{
  Assert(QuadratureSelector<dim>(name, order).get_points() == q.get_points(),
         ExcInternalError());
  deallog << name << ' ' << order << " ok" << std::endl;
}


int
main()
{
  initlog();

  check("gauss", 2, QGauss<1>(2));
  check("gauss", 2, QGauss<2>(2));
  check("gauss", 2, QGauss<3>(2));

  check("gauss", 2, QGauss<3>(2));
  check("gauss", 6, QGauss<3>(6));
  check("gauss", 10, QGauss<3>(10));

  check("weddle", 0, QWeddle<2>());

  check("gauss_lobatto", 2, QGaussLobatto<1>(2));
  check("gauss_lobatto", 3, QGaussLobatto<2>(3));
  check("gauss_lobatto", 4, QGaussLobatto<3>(4));

  check("gauss_chebyshev", 2, QGaussChebyshev<1>(2));
  check("gauss_chebyshev", 3, QGaussChebyshev<2>(3));
  check("gauss_chebyshev", 4, QGaussChebyshev<3>(4));

  check("gauss_radau_chebyshev", 2, QGaussRadauChebyshev<1>(2));
  check("gauss_radau_chebyshev", 3, QGaussRadauChebyshev<2>(3));
  check("gauss_radau_chebyshev", 4, QGaussRadauChebyshev<3>(4));

  check("gauss_lobatto_chebyshev", 2, QGaussLobattoChebyshev<1>(2));
  check("gauss_lobatto_chebyshev", 3, QGaussLobattoChebyshev<2>(3));
  check("gauss_lobatto_chebyshev", 4, QGaussLobattoChebyshev<3>(4));
}
