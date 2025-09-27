// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check for something in
// DerivativeApproximation::approximate_gradient. the original report
// stated it was a bug, but it was not (see the archives), but since
// the program is already there let's make use of it.
//
// this program is a modified version of one by
// Denis Danilov <danilovdenis@yandex.ru>,


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


class F : public Function<2>
{
public:
  virtual double
  value(const Point<2> &p, const unsigned int) const
  {
    double delta = 0.05;
    double x, y, r;
    x = p[0];
    y = p[1];
    r = std::sqrt(x * x + y * y);
    return 0.5 * (1 - std::tanh((r - 0.5) / (2 * M_SQRT2 * delta)));
  }
};


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  Triangulation<2> triangulation;
  FE_Q<2>          fe(2);
  DoFHandler<2>    dof_handler(triangulation);
  Vector<double>   phi_solution;
  Vector<float>    gradient_phi;
  float            gradient_phi_min, gradient_phi_max;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(5);

  dof_handler.distribute_dofs(fe);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();

  phi_solution.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, F(), phi_solution);

  gradient_phi.reinit(triangulation.n_active_cells());
  DerivativeApproximation::approximate_gradient(dof_handler,
                                                phi_solution,
                                                gradient_phi);

  gradient_phi_min = 1e30;
  gradient_phi_max = -1;

  cell = dof_handler.begin_active();
  for (unsigned int cell_no = 0; cell != endc; ++cell, ++cell_no)
    {
      if (gradient_phi(cell_no) < gradient_phi_min)
        gradient_phi_min = gradient_phi(cell_no);

      if (gradient_phi(cell_no) > gradient_phi_max)
        gradient_phi_max = gradient_phi(cell_no);
    }

  deallog << "gradient_phi_min: " << gradient_phi_min << std::endl;
  deallog << "gradient_phi_max: " << gradient_phi_max << std::endl;
}
