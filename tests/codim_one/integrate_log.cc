// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// integrates the function ln(x-a)*f(x), where f(x) is a power of x on the interval [a,b].
// dim=1 only.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files needed for the program

#include <deal.II/base/exceptions.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>


#include <fstream>
#include <string>

#include <math.h>

std::ofstream logfile("output");

double test(const unsigned int n,
            const unsigned int exponent,
            const double a,
            const double b)
{

  const QGaussLog<1> qlog(n);
  const QGauss<1> qgauss(7);

  Tensor<1,1> exp;
  exp[0]=exponent;
  Functions::Monomial<1> monomial(exp);

  Triangulation<1> tria;
  GridGenerator::hyper_cube(tria,a,b);

  FE_Q<1> feq(1);

  DoFHandler<1> dh(tria);

  FEValues<1> fe_values(feq, qlog, update_JxW_values | update_quadrature_points),
           fev_help(feq, qgauss, update_JxW_values | update_quadrature_points);


  dh.distribute_dofs(feq);

  double integral=0;

  DoFHandler<1>::active_cell_iterator cell=dh.begin_active(),
                                      endc=dh.end();

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit(cell);
      fev_help.reinit(cell);
      for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
        {
          integral+=
            fe_values.JxW(i)
            *
            monomial.value(fe_values.quadrature_point(i));
        }


      // The quadrature formula weights are chosen in order to
      // integrate f(x)*ln[(x-a)/(b-a)] (that is, the argument of the
      // logaritm goes linearly from 0 to 1 on the interval of
      // integration). So, the only thing one can specify is the
      // function f(x). Using logaritm properties, one can add the
      // integral of f(x)*ln(b-a) to this term to finally obtain the
      // integral of f(x)*ln(x-a).

      for (unsigned int i=0; i<fev_help.n_quadrature_points; ++i)
        {
          integral+=
            fev_help.JxW(i)
            *
            monomial.value(fev_help.quadrature_point(i))
            *
            log(b-a);

        }
    }

  return integral;

}

unsigned int
factorial(unsigned int a)
{
  if ( a == 0) return 1;

  if ( a-1 != 0 )
    a *= factorial(a-1);

  return a;
}

double
newton_binomial(unsigned int a,
                unsigned int b)
{
  double c;
  if (a >= b)
    c = factorial(a)
        / factorial(b)
        / factorial(a-b);
  else
    deallog<<"Error: in Newton binomial the first integer must be greater or equal to the second.";

  return c;
}

int main()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog<<std::fixed;
  deallog<<std::setprecision(5);
  const double a = 1,
               b = 5.;

  double exact_integral;

  deallog << "Calculation of the integral of ln(x-a)*f(x) on the interval [a,b] = ["
          << a << ", " << b << "]"<< std::endl;
  for (unsigned int j=0; j<13; ++j)
    {
      exact_integral=0;
      for (unsigned int k=0; k<=j; k++)
        exact_integral +=
          newton_binomial(j,k)
          *pow(a,static_cast<int>(j-k))
          *(
            pow(b-a,static_cast<int>(k)+1) / (k+1) * log(b-a)
            -
            pow(b-a,static_cast<int>(k)+1) / pow(static_cast<double>(k+1), 2)
          );
      deallog << "f(x) = x^" << j << std::endl;
      for (unsigned int i=1; i<13; ++i)
        deallog << " No. of points = " << i << "  "
                << "Error: "<< test(i,j,a,b)-exact_integral
                << std::endl;
      deallog<<std::endl;
    }
}
