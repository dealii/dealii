// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test is the same as matrix_vector_rt_01.cc but with hanging nodes.

#include <deal.II/dofs/dof_tools.h>

#include "../tests.h"

#include "matrix_vector_rt_common.h"


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(2);
  else
    tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active(), endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_RaviartThomasNodal<dim> fe(fe_degree - 1);
  DoFHandler<dim>            dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  deallog << "Using " << dof.get_fe().get_name() << std::endl;
  deallog << "Number of cells: " << dof.get_triangulation().n_active_cells()
          << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl
          << std::endl;
  do_test<dim, fe_degree, double>(dof, constraints, TestType::values_gradients);
}
