// ---------------------------------------------------------------------
// $Id: system_01.cc 31349 2013-10-20 19:07:06Z maier $
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// test get_name() (it is currently broken for FE_DGQArbitraryNodes)
/*
6264: An error occurred in line <833> of file </ssd/deal-trunk/deal.II/source/fe/fe_dgq.cc> in function
6264:     std::string dealii::FE_DGQArbitraryNodes<dim, spacedim>::get_name() const [with int dim = 2, int spacedim = 2, std::string = std::basic_string<char>]
6264: The violated condition was: 
6264:     index == n_points
6264: The name and call sequence of the exception was:
6264:     ExcMessage ("Could not decode support points in one coordinate direction.")
6264: Additional Information: 
6264: Could not decode support points in one coordinate direction.
6264: --------------------------------------------------------
*/
  
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


template <int dim>
void test(const FiniteElement<dim> &fe)
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
    QGauss<1> quadrature_g(5);
    FE_DGQArbitraryNodes<2> fe(quadrature_g);
    test(fe);
  }
  {
    QGaussLobatto<1> quadrature_gl(5);
    FE_DGQArbitraryNodes<2> fe(quadrature_gl);
    test(fe);
  }
  {
    QGaussLog<1> quadrature(3);
    FE_DGQArbitraryNodes<2> fe(quadrature);
    test(fe);
  }

  

  return 0;
}



