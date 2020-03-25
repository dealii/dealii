// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Test Legendre expansion in 2D and 3D for a function given using Legendre
// coefficients.
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <gsl/gsl_sf_legendre.h>

#include <iostream>

#include "../tests.h"


template <int dim>
class LegendreFunction : public Function<dim>
{
public:
  LegendreFunction(const Table<dim, double> &coefficients)
    : dealii::Function<dim>(1)
    , coefficients(coefficients)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const;

  const Table<dim, double> &
  get_coefficients() const
  {
    return coefficients;
  }

private:
  const Table<dim, double> coefficients;
};

// copy-paste from fe_series.cc
template <int dim>
double
Lh(const Point<dim> &x_q, const TableIndices<dim> &indices)
{
  double res = 1.0;
  for (unsigned int d = 0; d < dim; d++)
    {
      const double x = 2.0 * (x_q[d] - 0.5);
      Assert((x_q[d] <= 1.0) && (x_q[d] >= 0.),
             ExcMessage("x_q is not in [0,1]" + Utilities::to_string(x_q[d])));
      const int ind = indices[d];
      res *= sqrt(2.0) * gsl_sf_legendre_Pl(ind, x);
    }
  return res;
}

template <>
double
LegendreFunction<2>::value(const dealii::Point<2> &point,
                           const unsigned int) const
{
  double f = 0.0;

  for (unsigned int i = 0; i < coefficients.size(0); i++)
    for (unsigned int j = 0; j < coefficients.size(1); j++)
      f += Lh(point, TableIndices<2>(i, j)) * coefficients(i, j);

  return f;
}

template <>
double
LegendreFunction<3>::value(const dealii::Point<3> &point,
                           const unsigned int) const
{
  double f = 0.0;

  for (unsigned int i = 0; i < coefficients.size(0); i++)
    for (unsigned int j = 0; j < coefficients.size(1); j++)
      for (unsigned int k = 0; k < coefficients.size(2); k++)
        f += Lh(point, TableIndices<3>(i, j, k)) * coefficients(i, j, k);

  return f;
}

void
print(const Table<2, double> &coeff)
{
  for (unsigned int i = 0; i < coeff.size(0); i++)
    for (unsigned int j = 0; j < coeff.size(1); j++)
      deallog << coeff(i, j) << " ";
  deallog << std::endl;
}

void
print(const Table<3, double> &coeff)
{
  for (unsigned int i = 0; i < coeff.size(0); i++)
    for (unsigned int j = 0; j < coeff.size(1); j++)
      for (unsigned int k = 0; k < coeff.size(2); k++)
        deallog << coeff(i, j, k) << " ";
  deallog << std::endl;
}

void resize(Table<2, double> &coeff, const unsigned int N)
{
  coeff.reinit(N, N);
}

void resize(Table<3, double> &coeff, const unsigned int N)
{
  TableIndices<3> size;
  for (unsigned int d = 0; d < 3; d++)
    size[d] = N;
  coeff.reinit(size);
}



template <int dim>
void
test(const LegendreFunction<dim> &func, const unsigned int poly_degree)
{
  const unsigned int max_poly = poly_degree + 3;
  deallog << "-----------------------------------" << std::endl;
  deallog << dim << "d, p=" << poly_degree << ", max_p=" << max_poly
          << std::endl;
  deallog << "-----------------------------------" << std::endl;
  Triangulation<dim>    triangulation;
  hp::DoFHandler<dim>   dof_handler(triangulation);
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim>  quadrature_formula;

  // add some extra FEs in fe_collection
  for (unsigned int p = 1; p <= max_poly; p++)
    {
      fe_collection.push_back(FE_Q<dim>(p));
      quadrature_formula.push_back(QGauss<dim>(p + 1 + 5));
    }

  GridGenerator::hyper_cube(triangulation, 0.0, 1.0); // reference cell
  const unsigned int fe_index = poly_degree - 1;
  dof_handler.begin_active()->set_active_fe_index(fe_index);
  dof_handler.distribute_dofs(fe_collection);

  Vector<double> values(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, func, values);

  const unsigned int      N = poly_degree + 1;
  FESeries::Legendre<dim> legendre(N, fe_collection, quadrature_formula);

  const Table<dim, double> &coeff_in = func.get_coefficients();
  Table<dim, double>        coeff_out;
  resize(coeff_out, N);

  Vector<double> local_dof_values;

  typename hp::DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  {
    const unsigned int cell_n_dofs          = cell->get_fe().dofs_per_cell;
    const unsigned int cell_active_fe_index = cell->active_fe_index();

    local_dof_values.reinit(cell_n_dofs);
    cell->get_dof_values(values, local_dof_values);

    legendre.calculate(local_dof_values, cell_active_fe_index, coeff_out);
  }

  deallog << "calculated:" << std::endl;
  print(coeff_out);
  deallog << "exact:" << std::endl;
  print(coeff_in);

  dof_handler.clear();
}

int
main()
{
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  {
    const unsigned int dim      = 2;
    const unsigned int coeff_1d = 2;
    const unsigned int p        = 1;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; i++)
      for (unsigned int j = 0; j < coeff_1d; j++)
        coeff_in(i, j) = 1.0 + ind++;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);
  }

  {
    const unsigned int dim      = 2;
    const unsigned int coeff_1d = 3;
    const unsigned int p        = 2;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; i++)
      for (unsigned int j = 0; j < coeff_1d; j++)
        coeff_in(i, j) = 1.0 + ind++;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);
  }

  {
    const unsigned int dim      = 3;
    const unsigned int coeff_1d = 2;
    const unsigned int p        = 1;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; i++)
      for (unsigned int j = 0; j < coeff_1d; j++)
        for (unsigned int k = 0; k < coeff_1d; k++)
          coeff_in(i, j, k) = 1.0 + ind++;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);
  }

  {
    const unsigned int dim      = 3;
    const unsigned int coeff_1d = 3;
    const unsigned int p        = 2;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; i++)
      for (unsigned int j = 0; j < coeff_1d; j++)
        for (unsigned int k = 0; k < coeff_1d; k++)
          coeff_in(i, j, k) = 1.0 + ind++;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);
  }

  dealii::deallog << "Ok" << std::endl;
}
