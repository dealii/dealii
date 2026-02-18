// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


int
main()
{
  initlog();

  const int dim = 2;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe(2);

  GridGenerator::hyper_cube(triangulation);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  Vector<std::complex<double>> solution;
  solution.reinit(dof_handler.n_dofs());
  constraints.distribute(solution);

  deallog << "OK" << std::endl;
}
