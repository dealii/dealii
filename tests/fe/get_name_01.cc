// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test get_name()

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <string>

#include "../tests.h"


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
}


int
main()
{
  initlog();

  {
    FE_Q<2> fe(1);
    test(fe);
  }
  {
    FE_Q<2> fe(3);
    test(fe);
  }
  {
    FE_Q<2> fe(QIterated<1>(QTrapezoid<1>(), 3));
    test(fe);
  }
  {
    QGauss<1>               quadrature_g(5);
    FE_DGQArbitraryNodes<2> fe(quadrature_g);
    test(fe);
  }
  {
    QGaussLobatto<1>        quadrature_gl(5);
    FE_DGQArbitraryNodes<2> fe(quadrature_gl);
    test(fe);
  }
  {
    QGaussLog<1>            quadrature(3);
    FE_DGQArbitraryNodes<2> fe(quadrature);
    test(fe);
  }
  {
    QIterated<1>            quadrature(QTrapezoid<1>(), 3);
    FE_DGQArbitraryNodes<2> fe(quadrature);
    test(fe);
  }



  return 0;
}
