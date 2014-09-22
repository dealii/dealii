// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// Output function values and derivatives for vector valued
// polynomials returning Tensor<1,dim> for their values.

// classes: PolynomialsBDM, PolynomialsRaviartThomas

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_raviart_thomas.h>

#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

template<int dim, class POLY>
void check_point (const Point<dim> &x,
                  const POLY &p)
{
  const unsigned int n = p.n();
  std::vector<Tensor<1,dim> > values(n);
  std::vector<Tensor<2,dim> > gradients(n);
  std::vector<Tensor<3,dim> > seconds(0);

  p.compute(x, values, gradients, seconds);

  deallog << "Point " << x << std::endl;
  for (unsigned int i=0; i<n; ++i)
    {
      deallog << "p[" << i << "] value ";
      for (unsigned int d=0; d<dim; ++d)
        deallog << (int) (values[i][d]+.5) << ' ';
      deallog << " first ";
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          deallog << (int) (gradients[i][d1][d2]+.5) << ' ';
      deallog << std::endl;
    }
}


template<int dim>
void check_bdm ()
{
  Point<dim> x;

  PolynomialsBDM<dim> p1(1);
  PolynomialsBDM<dim> p2(2);
  PolynomialsBDM<dim> p3(3);
  PolynomialsBDM<dim> p4(4);

  x(0) = 2.;
  x(1) = 3.;
  if (dim>2)
    x(2) = 4;

  check_point(x, p1);
  check_point(x, p2);
  check_point(x, p3);
  check_point(x, p4);
}

template<int dim>
void check_rt ()
{
  Point<dim> x;

  PolynomialsRaviartThomas<dim> p0(0);
  PolynomialsRaviartThomas<dim> p1(1);
  PolynomialsRaviartThomas<dim> p2(2);
  PolynomialsRaviartThomas<dim> p3(3);

  x(0) = 2.;
  x(1) = 3.;
  if (dim>2)
    x(2) = 4;

  check_point(x, p0);
  check_point(x, p1);
  check_point(x, p2);
  check_point(x, p3);
}

int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("BDM");
  check_bdm<2>();
  check_bdm<3>();
  deallog.pop();
  deallog.push("RT");
  check_rt<2>();
  check_rt<3>();
  deallog.pop();
}
