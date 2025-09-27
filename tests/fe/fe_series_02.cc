// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Legendre expansion in 1D for a function given using Legendre functions.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx17/cmath.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


/**
 * A 1D function given by input legendre coefficients
 */
template <int dim>
class LegendreFunction : public Function<dim>
{
public:
  LegendreFunction(const std::vector<double> coefficients)
    : dealii::Function<dim>(1)
    , coefficients(coefficients)
  {}

  virtual double
  value(const dealii::Point<dim> &point,
        const unsigned int        component = 0) const;

  const std::vector<double> &
  get_coefficients() const
  {
    return coefficients;
  }

private:
  const std::vector<double> coefficients;
};

template <int dim>
double
LegendreFunction<dim>::value(const dealii::Point<dim> &point,
                             const unsigned int) const
{
  Assert(dim == 1, dealii::ExcNotImplemented());

  double f = 0.0;

  for (unsigned int l = 0; l < coefficients.size(); ++l)
    {
      const double m = 0.5;                // mid-point
      const double h = 0.5;                // half-length
      const double x = (point[0] - m) / h; // 1D only
      f += sqrt(1.0 / h) * std_cxx17::legendre(l, x) * coefficients[l];
    }

  return f;
}

template <int dim>
void
test(const LegendreFunction<dim> &func, const unsigned int poly_degree)
{
  Triangulation<dim>    triangulation;
  DoFHandler<dim>       dof_handler(triangulation);
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(dealii::FE_Q<dim>(poly_degree));

  hp::QCollection<dim> quadrature_formula;
  quadrature_formula.push_back(QGauss<dim>(poly_degree + 6));

  // reference cell:
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);

  dof_handler.distribute_dofs(fe_collection);

  Vector<double> values(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, func, values);

  const unsigned int              N = poly_degree + 1;
  const std::vector<unsigned int> n_coefficients_per_direction(
    fe_collection.size(), N);
  FESeries::Legendre<dim> legendre(n_coefficients_per_direction,
                                   fe_collection,
                                   quadrature_formula);

  const std::vector<double> &coeff_in = func.get_coefficients();
  Table<1, double>           coeff_out(N);

  Vector<double> local_dof_values;

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();

  {
    const unsigned int cell_n_dofs          = cell->get_fe().dofs_per_cell;
    const unsigned int cell_active_fe_index = cell->active_fe_index();

    local_dof_values.reinit(cell_n_dofs);
    cell->get_dof_values(values, local_dof_values);

    legendre.calculate(local_dof_values, cell_active_fe_index, coeff_out);
  }

  for (unsigned int i = 0; i < coeff_in.size(); ++i)
    deallog << coeff_in[i] << ' ';

  deallog << std::endl;

  for (unsigned int i = 0; i < N; ++i)
    deallog << coeff_out[i] << ' ';

  deallog << std::endl;

  dof_handler.clear();
}

int
main()
{
  const int dim = 1;

  initlog();

  {
    std::vector<double> coeff_in(2);
    coeff_in[0] = 1.0;
    coeff_in[1] = 2.0;
    LegendreFunction<dim> function(coeff_in);
    test(function, 1);
  }

  {
    std::vector<double> coeff_in(3);
    coeff_in[0] = 1.0;
    coeff_in[1] = 2.0;
    coeff_in[2] = 3.0;
    LegendreFunction<dim> function(coeff_in);
    test(function, 2);
  }

  {
    std::vector<double> coeff_in(4);
    coeff_in[0] = 1.0;
    coeff_in[1] = 2.0;
    coeff_in[2] = 3.0;
    coeff_in[3] = 4.0;
    LegendreFunction<dim> function(coeff_in);
    test(function, 3);
  }

  {
    std::vector<double> coeff_in(5);
    coeff_in[0] = 1.0;
    coeff_in[1] = 2.0;
    coeff_in[2] = 3.0;
    coeff_in[3] = 4.0;
    coeff_in[4] = 5.0;
    LegendreFunction<dim> function(coeff_in);
    test(function, 4);
  }

  {
    std::vector<double> coeff_in(6);
    coeff_in[0] = 1.0;
    coeff_in[1] = 2.0;
    coeff_in[2] = 3.0;
    coeff_in[3] = 4.0;
    coeff_in[4] = 5.0;
    coeff_in[5] = 6.0;
    LegendreFunction<dim> function(coeff_in);
    test(function, 5);
  }

  dealii::deallog << "Ok" << std::endl;
}
