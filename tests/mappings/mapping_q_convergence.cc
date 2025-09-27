// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <memory>

#include "../tests.h"

// Test that we achieve the correct rates of convergence with MappingQ when
// all cells are curved. Credit for this test case goes to Alexander
// Grayver. This test verifies that several issues with the combination of
// MappingQ and a ChartManifold have been fixed: more precisely, the computed
// errors for this program were, for a time, much lower for the 8.2 release
// than for the development branch. At the time of writing this (on 8.5pre)
// the results of this test are identical when compiled against either 8.2 or
// 8.5pre.
//
// In addition to checking the L2 errors of a projection, this test also
// verifies that there is no loss in convergence order.

// Like FESeries::linear_regression. Included here so that this test can also
// be run with older versions of deal.II.
double
regression_slope(const std::vector<double> &x, const std::vector<double> &y)
{
  FullMatrix<double> K(2, 2), invK(2, 2);
  Vector<double>     X(2), B(2);

  Assert(x.size() == y.size(),
         ExcMessage("x and y are expected to have the same size"));

  Assert(x.size() >= 2,
         dealii::ExcMessage(
           "at least two points are required for linear regression fit"));

  double sum_1 = 0.0, sum_x = 0.0, sum_x2 = 0.0, sum_y = 0.0, sum_xy = 0.0;

  for (unsigned int i = 0; i < x.size(); ++i)
    {
      sum_1 += 1.0;
      sum_x += x[i];
      sum_x2 += x[i] * x[i];
      sum_y += y[i];
      sum_xy += x[i] * y[i];
    }

  K(0, 0) = sum_1;
  K(0, 1) = sum_x;
  K(1, 0) = sum_x;
  K(1, 1) = sum_x2;

  B(0) = sum_y;
  B(1) = sum_xy;

  invK.invert(K);
  invK.vmult(X, B, false);

  return X(1);
}


double
zvalue(const double x, const double y)
{
  double xh = x * 5., yh = y * 5.;
  return (xh * exp(-xh * xh - yh * yh)) / 10.;
}

template <int dim>
class Geometry : public ChartManifold<dim>
{
public:
  virtual Point<dim>
  pull_back(const Point<dim> &space_point) const override;
  virtual Point<dim>
  push_forward(const Point<dim> &chart_point) const override;
  virtual std::unique_ptr<Manifold<dim>>
  clone() const override;
};

template <int dim>
Point<dim>
Geometry<dim>::pull_back(const Point<dim> &space_point) const
{
  const double d = space_point[dim - 1];
  const double z = zvalue(space_point[0], dim == 3 ? space_point[1] : 0);

  double d_hat = 0.;
  if ((d - z) <= 0)
    d_hat = (d - z) / (1. + z);
  else
    d_hat = (d - z) / (1. - z);

  Point<dim> p;
  for (unsigned i = 0; i < dim - 1; ++i)
    p[i] = space_point[i];
  p[dim - 1] = d_hat;

  return p;
}

template <int dim>
Point<dim>
Geometry<dim>::push_forward(const Point<dim> &chart_point) const
{
  const double d_hat = chart_point[dim - 1];
  const double z     = zvalue(chart_point[0], dim == 3 ? chart_point[1] : 0);

  double d = 0.;
  if (d_hat <= 0)
    d = d_hat + (d_hat + 1.) * z;
  else
    d = d_hat - (d_hat - 1.) * z;

  Point<dim> p;
  for (unsigned i = 0; i < dim - 1; ++i)
    p[i] = chart_point[i];
  p[dim - 1] = d;

  return p;
}

template <int dim>
std::unique_ptr<Manifold<dim>>
Geometry<dim>::clone() const
{
  return std::make_unique<Geometry<dim>>();
}


template <int dim>
class TranscendentalManufacturedSolution : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int /*component*/) const
  {
    return std::cos(p[0]) + 2.0 * std::sin(2 * p[1]);
  }
};



template <int dim>
void
create_tria(Triangulation<dim> &triangulation, const Geometry<dim> &geometry)
{
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  GridTools::transform(std::bind(&Geometry<dim>::push_forward,
                                 std::cref(geometry),
                                 std::placeholders::_1),
                       triangulation);

  triangulation.set_manifold(0, geometry);
  for (Triangulation<3>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    cell->set_all_manifold_ids(0);
}



template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  Geometry<dim>                           geometry;
  TranscendentalManufacturedSolution<dim> fe_function;
  AffineConstraints<double>               constraints;
  constraints.close();

  deallog << "FE degree: " << fe.degree << std::endl;
  for (unsigned mapping_p = 2; mapping_p < fe.degree + 3; ++mapping_p)
    {
      deallog << "mapping order: " << mapping_p << std::endl;
      Triangulation<dim> triangulation;
      create_tria(triangulation, geometry);
      DoFHandler<dim> dof_handler(triangulation);

      std::vector<double> log_refinements;
      std::vector<double> log_l2_errors;

      MappingQ<dim> mapping(mapping_p);
      for (unsigned int refinement_n = 1; refinement_n < 4; ++refinement_n)
        {
          triangulation.refine_global(1);
          dof_handler.clear();
          dof_handler.distribute_dofs(fe);

          Vector<double> v(dof_handler.n_dofs());
          VectorTools::project(mapping,
                               dof_handler,
                               constraints,
                               QGauss<dim>(fe.degree + (mapping_p + 2) / 2),
                               fe_function,
                               v);

          Vector<double> diff(triangulation.n_active_cells());
          VectorTools::integrate_difference(
            mapping,
            dof_handler,
            v,
            fe_function,
            diff,
            // superconvergence with QGauss(k + 1)
            // QGauss<dim>(fe.degree + 2),
            // normal convergence with QIterated
            QIterated<dim>(QGauss<1>(fe.degree + 1), 2),
            VectorTools::L2_norm);
          log_refinements.push_back(
            std::log10(std::pow(2.0, -double(refinement_n))));
          log_l2_errors.push_back(std::log10(diff.l2_norm()));
          if (log_refinements.size() > 1)
            {
              deallog << "current slope: "
                      << (*(log_l2_errors.end() - 1) -
                          *(log_l2_errors.end() - 2)) /
                           (*(log_refinements.end() - 1) -
                            *(log_refinements.end() - 2))
                      << std::endl;
            }
        }
      deallog << "Last number of DoFs: " << dof_handler.n_dofs() << std::endl;
      deallog << "Last L2 error: " << std::pow(10.0, log_l2_errors.back())
              << std::endl;
      deallog << "regression slope: "
              << regression_slope(log_refinements, log_l2_errors) << std::endl;
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(5);
  deallog.depth_console(0);

  const static unsigned dim = 3;

  for (unsigned p = 1; p < 4; ++p)
    {
      test<dim>(FE_Q<dim>(QGaussLobatto<1>(p + 1)));
    }
}
