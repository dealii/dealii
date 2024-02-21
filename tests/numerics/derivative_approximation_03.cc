// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Same as derivative_approximation_02, but without the mapping

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"


const double a = 4, b = 5;


template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    return sin(p[0] * a) * cos(p[1] * b);
  }
};

template <int dim>
void
exact_gradient(Point<dim> &p, Tensor<1, dim> &grad)
{
  double x = p[0], y = p[1];

  grad[0] = a * cos(a * x) * cos(b * y);
  grad[1] = -b * sin(a * x) * sin(b * y);
}

template <int dim>
void
exact_second(Point<dim> &p, Tensor<2, dim> &sec)
{
  double x = p[0], y = p[1];

  sec[0][0] = -a * a * sin(a * x) * cos(b * y);
  //  sec[0][1]=-a*b*cos(a*x)*sin(b*y);
  sec[1][0] = sec[0][1] = -a * b * cos(a * x) * sin(b * y);
  sec[1][1]             = -b * b * sin(a * x) * cos(b * y);
}

template <int dim>
void
exact_third(Point<dim> &p, Tensor<3, dim> &third)
{
  double x = p[0], y = p[1];
  // array of function and its derivatives
  double dx[4] = {sin(a * x),
                  a * cos(a * x),
                  -a * a * sin(a * x),
                  -a * a * a * cos(a * x)};
  double dy[4] = {cos(b * y),
                  -b * sin(b * y),
                  -b * b * cos(b * y),
                  b * b * b * sin(b * y)};

  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
        {
          int zeros = 0, ones = 0;
          switch (i)
            {
              case 0:
                ++zeros;
                break;
              case 1:
                ++ones;
                break;
              default:
                third[i][j][k] = 0;
                continue;
            }
          switch (j)
            {
              case 0:
                ++zeros;
                break;
              case 1:
                ++ones;
                break;
              default:
                third[i][j][k] = 0;
                continue;
            }
          switch (k)
            {
              case 0:
                ++zeros;
                break;
              case 1:
                ++ones;
                break;
              default:
                third[i][j][k] = 0;
                continue;
            }
          third[i][j][k] = dx[zeros] * dy[ones];
        }
}


template <int dim>
void
derivatives()
{
  MyFunction<dim>    function;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);
  FE_DGQ<dim>     fe(2);
  DoFHandler<dim> dof_handler(tria);
  Vector<double>  solution;
  QMidpoint<dim>  q_midpoint;
  FEValues<dim>   fe_values(fe, q_midpoint, update_quadrature_points);

  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, function, solution);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  for (; cell != endc; ++cell)
    {
      // get derivative approximations
      Tensor<1, dim> grad, ex_grad;
      Tensor<2, dim> second, ex_second;
      Tensor<3, dim> third, ex_third;

      DerivativeApproximation::approximate_derivative_tensor(
        dof_handler, solution, cell, grad, 0);
      double normgrad = DerivativeApproximation::derivative_norm(grad);

      DerivativeApproximation::approximate_derivative_tensor(
        dof_handler, solution, cell, second, 0);
      double normsecond = DerivativeApproximation::derivative_norm(second);

      DerivativeApproximation::approximate_derivative_tensor(
        dof_handler, solution, cell, third, 0);
      double normthird = DerivativeApproximation::derivative_norm(third);

      // symmetry of second derivative
      Assert(second == transpose(second),
             ExcMessage("Second derivative is not symmetric"));

      // symmetry of third derivative note,
      // that this is only part of the truth,
      // we would have to test more here to be
      // really sure, but this should be enough
      for (unsigned int i = 0; i < dim; ++i)
        Assert(third[i] == transpose(third[i]),
               ExcMessage("Third derivative is not symmetric"));

      // get exact derivatives
      fe_values.reinit(cell);
      std::vector<Point<dim>> q_point;
      q_point = fe_values.get_quadrature_points();

      exact_gradient(q_point[0], ex_grad);
      exact_second(q_point[0], ex_second);
      exact_third(q_point[0], ex_third);

      double ex_normgrad = DerivativeApproximation::derivative_norm(ex_grad);
      double ex_normsecond =
        DerivativeApproximation::derivative_norm(ex_second);
      double ex_normthird = DerivativeApproximation::derivative_norm(ex_third);

      // output of all values for comparison
      deallog << "cell " << cell << std::endl;
      deallog << "approx. gradient: " << grad << " , norm: " << normgrad
              << std::endl
              << "  exact gradient: " << ex_grad << " , norm: " << ex_normgrad
              << std::endl;
      deallog << "approx. second derivative (hessian): " << second
              << " , norm: " << normsecond << std::endl
              << "  exact second derivative (hessian): " << ex_second
              << " , norm: " << ex_normsecond << std::endl;
      deallog << "approx. third derivative: " << third
              << " , norm: " << normthird << std::endl
              << "  exact third derivative: " << ex_third
              << " , norm: " << ex_normthird << std::endl
              << std::endl;
    }
}


int
main()
{
  initlog();

  deallog << "------------ 2D ------------" << std::endl << std::endl;
  derivatives<2>();
  deallog << std::endl
          << "------------ 3D ------------" << std::endl
          << std::endl;
  derivatives<3>();
}
