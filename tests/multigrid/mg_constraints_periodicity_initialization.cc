// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors

// Check the default periodic constraints, the explicit opt-out, and resetting
// an already initialized object independently of the rotated transfer tests.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"

int
main()
{
  initlog();
  Triangulation<2> tria(Triangulation<2>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria, 0., 1., true);
  std::vector<GridTools::PeriodicFacePair<Triangulation<2>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(tria, 0, 1, 0, periodic_faces);
  tria.add_periodicity(periodic_faces);
  tria.refine_global(1);
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center()[1] > 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<2> dofs(tria);
  dofs.distribute_dofs(FE_Q<2>(1));
  dofs.distribute_mg_dofs();
  MGLevelObject<IndexSet> relevant(0, tria.n_global_levels() - 1);
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    relevant[level] =
      DoFTools::extract_locally_relevant_level_dofs(dofs, level);

  MGConstrainedDoFs defaults, with_relevant, explicit_option;
  defaults.initialize(dofs);
  with_relevant.initialize(dofs, relevant);
  explicit_option.initialize(dofs, relevant, false);

  const auto check = [&](const MGConstrainedDoFs &mg, const bool periodic) {
    for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
      {
        AssertThrow(mg.get_refinement_edge_indices(level) ==
                      defaults.get_refinement_edge_indices(level),
                    ExcInternalError());
        const auto &expected = defaults.get_level_constraints(level);
        const auto &actual   = mg.get_level_constraints(level);
        AssertThrow(expected.n_constraints() > 0, ExcInternalError());
        if (!periodic)
          AssertThrow(actual.n_constraints() == 0, ExcInternalError());
        else
          for (const auto i : relevant[level])
            {
              AssertThrow(actual.is_constrained(i) ==
                            expected.is_constrained(i),
                          ExcInternalError());
              if (expected.is_constrained(i))
                {
                  AssertThrow(*actual.get_constraint_entries(i) ==
                                *expected.get_constraint_entries(i),
                              ExcInternalError());
                  AssertThrow(actual.get_inhomogeneity(i) ==
                                expected.get_inhomogeneity(i),
                              ExcInternalError());
                }
            }
      }
  };

  check(with_relevant, true);
  check(explicit_option, false);
  explicit_option.initialize(dofs, relevant, true);
  check(explicit_option, true);
  // Reinitialization must remove the previous identity constraints.
  explicit_option.initialize(dofs, relevant, false);
  check(explicit_option, false);
  explicit_option.initialize(dofs);
  check(explicit_option, true);
  deallog << "OK" << std::endl;
}
