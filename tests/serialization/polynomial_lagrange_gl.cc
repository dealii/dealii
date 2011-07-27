//----------------------------  polynomial_lagrange_gl.cc  ---------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  polynomial_lagrange_gl.cc  ---------------

// check serialization for Lagrange polynomial based on Gauss-Lobatto support
// points

#include "serialization.h"
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  unsigned int n1 = 3;
  unsigned int support_point1 = 1;
  QGaussLobatto<1> quad1 (n1+1);
  Polynomials::Polynomial<double> p1 = Polynomials::
    generate_complete_Lagrange_basis (quad1.get_points())[support_point1];

  unsigned int n2 = 4;
  unsigned int support_point2 = 2;
  QGaussLobatto<1> quad2 (n2+1);
  Polynomials::Polynomial<double> p2 = Polynomials::
    generate_complete_Lagrange_basis (quad2.get_points())[support_point2];

  verify (p1, p2);
}


int main ()
{
  std::ofstream logfile("polynomial_lagrange_gl/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
