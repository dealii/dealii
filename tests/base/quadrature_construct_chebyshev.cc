// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// check that the QGaussChebyshev, QGaussRadauChebyshev and QGaussLobattoChebyshev,
// can be constructed in all dimensions. Previously, this failed since the base class
// constructor used required the one-dimensional quadrature formula to integrate
// constants exactly. This is not true for the classes considered here.


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
