// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// Test Legendre expansion in 1D for quadratic function coming from FE.
// Also test that our interpretation of GSL function is correct to have
// an orthogonal basis.

// MWE in Maxima
/**************************************
hj : 1/2;
mj : 1/2;
define(f(x), (1.81735e-05*(1.0-x)*(0.5-x)*2 + 0.000901649*x*(x-0.5)*2 + 1.35059e-05*x*(1.0-x)*4.0));
load("orthopoly");
orthopoly_returns_intervals : false;
plot2d([legendre_p(0,x), legendre_p(1,x),legendre_p(2,x)], [x,-1,1])$
define(Lh(n,h,m,x), sqrt(1/h)*legendre_p(n,(x-m)/h));
define(C(n),integrate(f(x)*Lh(n,hj,mj,x),x,0,1)*(n+1/2));
bfloat(C(0)), nouns;
bfloat(C(1)), nouns;
bfloat(C(2)), nouns;
bfloat(C(3)), nouns;
 **************************************/

#include "../tests.h"
#include <iostream>
#include <fstream>

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/numerics/vector_tools.h>

#include <gsl/gsl_sf_legendre.h>

using namespace dealii;

template<int dim>
class LegendreFunction : public Function<dim>
{
public:
  LegendreFunction()
    :
    Function<dim>(1)
  {}

  virtual double value(const dealii::Point<dim> &point,
                       const unsigned int component = 0 ) const;
};

template<int dim>
double LegendreFunction<dim>::value(const Point<dim> &point,
                                    const unsigned int ) const
{
  Assert(dim==1,
         dealii::ExcNotImplemented());

  const double &x = point[0];
  return 1.81735e-05*(1.0-x)*(0.5-x)*2 + 0.000901649*x*(x-0.5)*2 + 1.35059e-05*x*(1.0-x)*4.0;
}


template<int dim>
void test(const LegendreFunction<dim> &func,
          const unsigned int poly_degree)
{
  Triangulation<dim> triangulation;
  hp::DoFHandler<dim> dof_handler(triangulation);
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quadrature_formula;

  for (unsigned int p = poly_degree; p<= poly_degree+3; p++)
    {
      fe_collection.push_back(dealii::FE_Q<dim>(p));
      quadrature_formula.push_back(dealii::QGauss<dim>(p+1+5));
    }

  // reference cell
  GridGenerator::hyper_cube (triangulation,0.0,1.0);

  dof_handler.distribute_dofs (fe_collection);

  Vector<double> values(dof_handler.n_dofs());

  VectorTools::interpolate (dof_handler,func,values);
  const unsigned int N = 4;
  FESeries::Legendre<dim> legendre(N,fe_collection,quadrature_formula);

  Table<1,double>     coeff_out(N);
  Vector<double> local_dof_values;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active();
  {
    const unsigned int cell_n_dofs = cell->get_fe().dofs_per_cell;
    const unsigned int cell_active_fe_index = cell->active_fe_index();

    local_dof_values.reinit (cell_n_dofs);
    cell->get_dof_values (values, local_dof_values);

    legendre.calculate(local_dof_values,
                       cell_active_fe_index,
                       coeff_out);

    deallog << "local dofs:";
    for (unsigned int i = 0; i < cell_n_dofs; i++)
      dealii::deallog << " " <<local_dof_values[i];

    dealii::deallog << std::endl;
  }

  deallog << "calculated:"<<std::endl;
  for (unsigned int i = 0; i < N; i++)
    deallog << coeff_out[i] << std::endl;

  std::vector<double> coeff_exp(3);
  // coeff calculated in maxima (see MWE above):
  coeff_exp[0] = 1.147688635236788e-4;
  coeff_exp[1] = 3.123557585310879e-4;
  coeff_exp[2] = 2.104375000953028e-4;
  deallog << "exact:"<<std::endl;
  for (unsigned int i = 0; i < coeff_exp.size(); i++)
    deallog << coeff_exp[i] << std::endl;

  dof_handler.clear();
}


/**
 * Small test to first output Legendre coefficients from GSL at -1,0,1
 * and then check that they are orthonormal
 */
void test_legendre_orthonormal(const unsigned int N)
{
  const unsigned int dim = 1;
  deallog << "Pl @ -1;0;1"<<std::endl;
  for (unsigned int l = 0; l < N; l++)
    {
      deallog << "l="<<l<<": ";
      for (double x = -1.0; x <=1.0; x+=1.0)
        deallog<< gsl_sf_legendre_Pl (l, x) << " ";

      deallog<<std::endl;
    }

  QGauss<dim> quadrature (8);
  deallog <<"orthogonality: " << std::endl;
  for (int k1 = 0; k1 < N; k1++)
    for (int k2 = 0; k2 < N; k2++)
      {
        double ortho = 0;
        for (unsigned int q=0; q<quadrature.size(); ++q)
          {
            const Point<dim> &x_q = quadrature.point(q);
            const double       m = 0.5; // mid-point
            const double       h = 0.5; // half-length
            const double       x = (x_q[0]-m)/h; // 1D only
            Assert (std::fabs(x) < 1.0,
                    dealii::ExcInternalError());
            const double L1 = std::sqrt(1.0/h) * gsl_sf_legendre_Pl (k1, x);
            const double L2 = std::sqrt(1.0/h) * gsl_sf_legendre_Pl (k2, x);
            ortho += L1 * L2 * quadrature.weight(q);
          }
        ortho *=(1.0+k1+k2)/2.0;

        deallog << "("<<k1<<","<<k2<<") = " <<ortho<<std::endl;
      }
  deallog << std::endl;
}

int main ()
{
  const int dim = 1;

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1e-8);

  test_legendre_orthonormal(3);
  LegendreFunction<dim> function;
  test(function,2);

  dealii::deallog << "Ok"<<std::endl;
}
