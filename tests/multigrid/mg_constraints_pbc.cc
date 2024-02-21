/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Test that DoFTools::make_periodicity_constraints() and MGConstrainedDoFs
 * creates the same constrains for periodic constraints in the case of
 * globally refined meshes.
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int dim           = 2;
  const unsigned int n_refinements = 2;

  Triangulation<dim> tria(
    Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria, 0, 1, true);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;

  GridTools::collect_periodic_faces(tria, 0, 1, 0, periodic_faces);

  tria.add_periodicity(periodic_faces);
  tria.refine_global(n_refinements);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(1));
  dof_handler.distribute_mg_dofs();

  MGConstrainedDoFs mg_constraints;
  mg_constraints.initialize(dof_handler);

  const auto &af_level = mg_constraints.get_level_constraints(n_refinements);

  AffineConstraints<double> af;
  DoFTools::make_periodicity_constraints<dim, dim, double>(
    dof_handler, 0, 1, 0, af);

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      Assert(af_level.is_constrained(i) == af.is_constrained(i),
             ExcInternalError());

      if (af_level.is_constrained(i) == false)
        continue;

      Assert(*af_level.get_constraint_entries(i) ==
               *af.get_constraint_entries(i),
             ExcInternalError());
    }

#if 0
  af_level.print(std::cout);
  std::cout << std::endl;
  af.print(std::cout);
  std::cout << std::endl;
#endif

  deallog << "OK!" << std::endl;
}
