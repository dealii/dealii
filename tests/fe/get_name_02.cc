// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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



// test get_name()

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature_lib.h>

#include <fstream>
#include <string>


template <int dim, int spacedim>
void test(const FiniteElement<dim,spacedim> &fe)
{
  deallog << fe.get_name() << std::endl;
}


int
main()
{
  initlog();

  {
    FE_Q<2,3> fe(1);
    test(fe);
  }
  {
    FE_Q<2,3> fe(QGaussLobatto<1>(4));
    test(fe);
  }
  {
    QGauss<1> quadrature_g(5);
    FE_DGQArbitraryNodes<2,3> fe(quadrature_g);
    test(fe);
  }
  {
    QGaussLobatto<1> quadrature_gl(5);
    FE_DGQArbitraryNodes<2,3> fe(quadrature_gl);
    test(fe);
  }
  {
    QGaussLog<1> quadrature(3);
    FE_DGQArbitraryNodes<2,3> fe(quadrature);
    test(fe);
  }
  {
    QIterated<1> quadrature(QTrapez<1>(), 3);
    FE_DGQArbitraryNodes<2,3> fe(quadrature);
    test(fe);
  }



  return 0;
}
