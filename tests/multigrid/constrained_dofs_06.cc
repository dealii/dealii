// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check make_no_normal_flux_constraints function for a box mesh

#include <deal.II/base/exceptions.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <vector>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(triangulation, 0., 1., true);
  triangulation.refine_global(4 - dim);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(1), dim));
  dof_handler.distribute_mg_dofs();

  std::set<types::boundary_id> no_flux_boundary = {0, 1, 2};

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(dof_handler);
  for (auto bid : no_flux_boundary)
    mg_constrained_dofs.make_no_normal_flux_constraints(dof_handler, bid, 0);

  for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
    {
      deallog << "Level " << level << ": " << std::flush;
      mg_constrained_dofs.get_boundary_indices(level).print(
        deallog.get_file_stream());
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  deallog << "2D" << std::endl;
  test<2>();
  deallog << std::endl << "3D" << std::endl;
  test<3>();
  return 0;
}
