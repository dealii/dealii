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



// check that the QGaussChebyshev, QGaussRadauChebyshev and
// QGaussLobattoChebyshev, can be constructed in all dimensions. Previously,
// this failed since the base class constructor used required the
// one-dimensional quadrature formula to integrate constants exactly. This is
// not true for the classes considered here.


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

template <int dim>
void
construct_quadrature()
{
  QGaussChebyshev<dim> q_1(2);
  deallog << "QGaussChebyshev<" << dim << ">: OK" << std::endl;
  QGaussRadauChebyshev<dim> q_2(2);
  deallog << "QGaussRadauChebyshev<" << dim << ">: OK" << std::endl;
  QGaussLobattoChebyshev<dim> q_3(2);
  deallog << "QGaussLobattoChebyshev<" << dim << ">: OK" << std::endl;
}

int
main()
{
  initlog();

  construct_quadrature<1>();
  construct_quadrature<2>();
  construct_quadrature<3>();
}
