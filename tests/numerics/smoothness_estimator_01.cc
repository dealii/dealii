// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// essentially similar to fe/fe_series_05.cc but test smoothness estimation.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx17/cmath.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/vector_tools.h>

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
  for (unsigned int d = 0; d < dim; ++d)
    {
      const double x = 2.0 * (x_q[d] - 0.5);
      Assert((x_q[d] <= 1.0) && (x_q[d] >= 0.),
             ExcMessage("x_q is not in [0,1]" + Utilities::to_string(x_q[d])));
      const unsigned int ind = indices[d];
      res *= sqrt(2.0) * std_cxx17::legendre(ind, x);
    }
  return res;
}



template <>
double
LegendreFunction<2>::value(const dealii::Point<2> &point,
                           const unsigned int) const
{
  double f = 0.0;

  for (unsigned int i = 0; i < coefficients.size(0); ++i)
    for (unsigned int j = 0; j < coefficients.size(1); ++j)
      f += Lh(point, TableIndices<2>(i, j)) * coefficients(i, j);

  return f;
}

template <>
double
LegendreFunction<3>::value(const dealii::Point<3> &point,
                           const unsigned int) const
{
  double f = 0.0;

  for (unsigned int i = 0; i < coefficients.size(0); ++i)
    for (unsigned int j = 0; j < coefficients.size(1); ++j)
      for (unsigned int k = 0; k < coefficients.size(2); ++k)
        f += Lh(point, TableIndices<3>(i, j, k)) * coefficients(i, j, k);

  return f;
}



void
compare(const Table<2, double> &coeff1, const Table<2, double> &coeff2)
{
  double linf = 0.;
  for (unsigned int i = 0; i < coeff1.size(0); ++i)
    for (unsigned int j = 0; j < coeff1.size(1); ++j)
      linf = std::max(linf, std::abs(coeff1(i, j) - coeff2(i, j)));

  deallog << "Linf norm in exact and calculate Legendre coefficients:"
          << std::endl
          << linf << std::endl;
}

void
compare(const Table<3, double> &coeff1, const Table<3, double> &coeff2)
{
  double linf = 0.;
  for (unsigned int i = 0; i < coeff1.size(0); ++i)
    for (unsigned int j = 0; j < coeff1.size(1); ++j)
      for (unsigned int k = 0; k < coeff1.size(2); ++k)
        linf = std::max(linf, std::abs(coeff1(i, j, k) - coeff2(i, j, k)));

  deallog << "Linf norm in exact and calculate Legendre coefficients:"
          << std::endl
          << linf << std::endl;
}



void
resize(Table<2, double> &coeff, const unsigned int N)
{
  coeff.reinit(N, N);
}

void
resize(Table<3, double> &coeff, const unsigned int N)
{
  TableIndices<3> size;
  for (unsigned int d = 0; d < 3; ++d)
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

  // add some extra FEs in fe_collection
  hp::FECollection<dim> fe_collection;
  for (unsigned int p = 1; p <= max_poly; ++p)
    fe_collection.push_back(FE_Q<dim>(p));

  FESeries::Legendre<dim> legendre =
    SmoothnessEstimator::Legendre::default_fe_series(fe_collection);

  const unsigned int fe_index = poly_degree - 1;
  const unsigned int n_modes =
    legendre.get_n_coefficients_per_direction(fe_index);

  // custom predicate:
  // p-ref for linear elements and use j=1,...,pe otherwise.
  ComponentMask coefficients_predicate(n_modes, true);
  coefficients_predicate.set(0, false);

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0); // reference cell

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.begin_active()->set_active_fe_index(fe_index);
  dof_handler.distribute_dofs(fe_collection);

  Vector<double> values(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, func, values);

  const Table<dim, double> &coeff_in = func.get_coefficients();
  Table<dim, double>        coeff_out;
  resize(coeff_out, n_modes);

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

  compare(coeff_in, coeff_out);

  // finally test smoothness estimator:
  Vector<float> smoothness(1);
  SmoothnessEstimator::Legendre::coefficient_decay_per_direction(
    legendre,
    dof_handler,
    values,
    smoothness,
    coefficients_predicate,
    /*smallest_abs_coefficient=*/1e-10,
    /*only_flagged_cells=*/false);

  deallog << "smoothness:" << std::endl << smoothness[0] << std::endl;

  dof_handler.clear();
}



int
main()
{
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  // for linear elements we expect p-refinement by convention
  {
    const unsigned int dim      = 2;
    const unsigned int coeff_1d = 2;
    const unsigned int p        = 1;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; ++i)
      for (unsigned int j = 0; j < coeff_1d; ++j)
        coeff_in(i, j) = 1.0 + ind++;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);
    deallog << "expected smoothness:" << std::endl
            << std::numeric_limits<float>::infinity() << std::endl;
  }

  // for quadratic we can already assign exponential decay:   a_i = C exp ( -k
  // i) set one with different k's
  {
    const double k1 = 1.;
    const double k2 = 2.;

    const unsigned int dim      = 2;
    const unsigned int coeff_1d = 3;
    const unsigned int p        = 2;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; ++i)
      coeff_in(i, 0) = exp(-k1 * i);

    for (unsigned int i = 0; i < coeff_1d; ++i)
      coeff_in(0, i) = exp(-k2 * i);

    // make sure predicate skips 0-th:
    coeff_in(0, 0) = 12345;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);

    deallog << "expected smoothness:" << std::endl
            << std::min(k1, k2) << std::endl;
  }

  // linear elements in 3D (expect zero output)
  {
    const unsigned int dim      = 3;
    const unsigned int coeff_1d = 2;
    const unsigned int p        = 1;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d, coeff_1d);
    unsigned int       ind = 0;
    for (unsigned int i = 0; i < coeff_1d; ++i)
      for (unsigned int j = 0; j < coeff_1d; ++j)
        for (unsigned int k = 0; k < coeff_1d; ++k)
          coeff_in(i, j, k) = 1.0 + ind++;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);
    deallog << "expected smoothness:" << std::endl
            << std::numeric_limits<float>::infinity() << std::endl;
  }

  // cubic in 3D
  {
    const double       k1       = 2.;
    const double       k2       = 3.;
    const double       k3       = 4.;
    const unsigned int dim      = 3;
    const unsigned int coeff_1d = 4;
    const unsigned int p        = 3;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d, coeff_1d);
    for (unsigned int i = 0; i < coeff_1d; ++i)
      coeff_in(i, 0, 0) = exp(-k1 * i);

    for (unsigned int j = 0; j < coeff_1d; ++j)
      coeff_in(0, j, 0) = exp(-k2 * j);

    for (unsigned int k = 0; k < coeff_1d; ++k)
      coeff_in(0, 0, k) = exp(-k3 * k);

    // make sure predicate skips 0-th:
    coeff_in(0, 0, 0) = 12345;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);

    deallog << "expected smoothness:" << std::endl
            << std::min(k1, std::min(k2, k3)) << std::endl;
  }


  // 4-th order in 3D but with some coefficients being zero
  {
    const double       k1       = 2.;
    const double       k2       = k1 + 1.;
    const unsigned int dim      = 3;
    const unsigned int coeff_1d = 5;
    const unsigned int p        = 4;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d, coeff_1d);
    // all non-zero:
    for (unsigned int i = 0; i < coeff_1d; ++i)
      coeff_in(i, 0, 0) = exp(-k2 * i);

    // some non-zero (2nd and 4th), the slowest decay will be from this
    // direction
    for (unsigned int j = 2; j < coeff_1d; j = j + 2)
      coeff_in(0, j, 0) = exp(-k1 * j);

    // all but one zero:
    for (unsigned int k = 3; k < coeff_1d; k = k + 10)
      coeff_in(0, 0, k) = exp(-k2 * k);

    // make sure predicate skips 0-th:
    coeff_in(0, 0, 0) = 12345;

    LegendreFunction<dim> function(coeff_in);
    test(function, p);

    deallog << "expected smoothness:" << std::endl << k1 << std::endl;
  }

  // cubic in 3D (zero)
  {
    const unsigned int dim      = 3;
    const unsigned int coeff_1d = 4;
    const unsigned int p        = 3;
    Table<dim, double> coeff_in(coeff_1d, coeff_1d, coeff_1d);

    LegendreFunction<dim> function(coeff_in);
    test(function, p);

    deallog << "expected smoothness:" << std::endl
            << std::numeric_limits<float>::infinity() << std::endl;
  }

  dealii::deallog << "Ok" << std::endl;
}
